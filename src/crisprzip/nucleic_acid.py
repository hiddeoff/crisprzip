"""Represents nucleic acid hybrids, either by mismatch positions or sequences."""

import os

import json
import random
from pathlib import Path
from typing import Union, List, Tuple, Callable

import numpy as np
from functools import lru_cache
from numpy.random import Generator, default_rng
from numpy.typing import ArrayLike


def format_point_mutations(protospacer: str,
                           target_sequence: str) -> List[str]:
    """List the point mutations between ``target_sequence`` and ``protospacer``."""

    if len(protospacer) != len(target_sequence):
        raise ValueError("Protospacer and target should be equally long.")

    mm_list = []
    for i, b in enumerate(protospacer):
        if target_sequence[i] == b:
            continue

        mm_list += [f"{b}{i+1:02d}{target_sequence[i]}"]

    return mm_list


class MismatchPattern:
    """Positions of the mismatched bases bases in a target sequence.

    Attributes
    ----------
    pattern : `numpy.ndarray`
        Array with True indicating mismatched basepairs
    length : `int`
        Guide length
    mm_num : `int`
        Number of mismatches in the array
    is_on_target : `bool`
        Indicates whether the array is the on-target array

    Notes
    -----
    Assumes a 3'-to-5' DNA direction. (CRISPR-Cas9 directionality).
    """

    def __init__(self, array: np.typing.ArrayLike):
        array = np.array(array)
        if array.ndim != 1:
            raise ValueError('Array should be 1-dimensional')
        if not (np.all((array == 0) | (array == 1)) or
                np.all((array is False) | (array is True)) or
                np.all((np.isclose(array, 0.0)) | (np.isclose(array, 0.0)))):
            raise ValueError('Array should only contain 0 and 1 values')

        self.pattern = np.asarray(array, dtype='bool')
        self.length = self.pattern.size
        self.mm_num = int(np.sum(self.pattern))
        self.is_on_target = self.mm_num == 0

    def __repr__(self):
        return "".join(["1" if mm else "0" for mm in self.pattern])

    def __str__(self):
        return self.__repr__()

    @classmethod
    def from_string(cls, mm_array_string):
        return cls(np.array(list(mm_array_string), dtype='int'))

    @classmethod
    def from_mm_pos(cls, guide_length: int, mm_pos_list: list = None,
                    zero_based_index=False):
        """Alternative constructor. Uses 1-based indexing by default."""
        array = np.zeros(guide_length)

        if mm_pos_list is None:
            mm_pos_list = []

        if not zero_based_index:
            mm_pos_list = [x - 1 for x in mm_pos_list]

        if mm_pos_list is not None:
            array[mm_pos_list] = 1
        return cls(array)

    @classmethod
    def from_target_sequence(cls, protospacer: str,
                             target_sequence: str) -> 'MismatchPattern':
        """Alternative constructor"""
        pmut_list = format_point_mutations(protospacer, target_sequence)
        return cls.from_mm_pos(
            len(protospacer),
            [int(pmut[1:3]) for pmut in pmut_list]
        )

    @classmethod
    def make_random(cls, guide_length: int, mm_num: int,
                    rng: Union[int, Generator] = None):
        if type(rng) is int or rng is None:
            rng = default_rng(rng)
        target = np.zeros(guide_length)
        mm_pos = rng.choice(range(20), size=mm_num, replace=False).tolist()
        target[mm_pos] = 1
        return cls(target)

    def get_mm_pos(self):
        return [i for i, mm in enumerate(self.pattern) if mm]


def make_hybr_energy_func(protospacer: str,
                          weight: Union[float, Tuple[float, float]] = None) \
        -> Callable:
    """Make a hybridization energy function."""

    def get_hybr_energy_fixed_protospacer(offtarget_seq: str):
        return get_hybridization_energy(protospacer, offtarget_seq,
                                        weight=weight)

    return get_hybr_energy_fixed_protospacer


def get_hybridization_energy(protospacer: str,
                             offtarget_seq: str = None,
                             mutations: str = '',
                             weight: Union[float,
                                           Tuple[float, float]] = None) \
        -> np.ndarray:
    """Calculate the free energy cost of R-loop formation.

    Parameters
    ----------
    protospacer : `str`
        Full sequence of the protospacer/on-target: 5'-20nt-PAM-3',
        possibly preceded by the upstream sequence.
    offtarget_seq : `str`
        Full sequence of the (off-)target: 5'-20nt-PAM-3', possibly
        preceded by the upstream sequence. If provided, overrules
        the ``mutations`` parameter.
    mutations : `str`
        Mismatch desciptors (in the form "A02T") describing how the
        target deviates from the protospacer. Multiple mismatches
        should be space-separated. Is empty by default, indicating
        no mismatches (=on-target hybridization energy).
    weight : `float` or `tuple`[`float`], optional
        Optional weighing of the dna opening energy and rna duplex energy.
        If `None` (default), no weighing is applied. If `float`, both DNA and
        RNA energies are multiplied by the weight parameter. If `tuple``
        of two `float`s, the first value is used as a multiplier for the
        DNA opening energy, and the second is used as a multiplier for the
        RNA-DNA hybridization energy.


    Returns
    -------
    hybridization_energy : `numpy.ndarray`
        Free energies required to create an R-loop.
    """

    # Handling 'mutation' argument
    if offtarget_seq is None:
        hybrid = GuideTargetHybrid.from_cas9_protospacer(protospacer,
                                                         mutations)
        # Recursive calling to include in caching
        offtarget_seq = hybrid.target.seq2 + protospacer[-3:]
        if hybrid.target.upstream_nt is not None:
            offtarget_seq = hybrid.target.upstream_nt + offtarget_seq
        return get_hybridization_energy(
            protospacer=protospacer,
            offtarget_seq=offtarget_seq,
            weight=weight
        )

    # do calculations
    out = get_na_energies_cached(protospacer, offtarget_seq)
    dna_energy = np.array(out[0])
    rna_energy = np.array(out[1])
    if weight is None:
        return dna_energy + rna_energy
    elif isinstance(weight, (float, int, np.floating, np.integer)):
        return weight * (dna_energy + rna_energy)
    elif type(weight) is tuple:
        return (weight[0] * dna_energy +
                weight[1] * rna_energy)


@lru_cache
def get_na_energies_cached(protospacer: str, offtarget_seq: str = None) -> \
        Tuple[Tuple[float, ], Tuple[float, ]]:
    """Calculate the DNA and RNA contributoins to the R-loop cost (with caching).

    Parameters
    ----------
    protospacer : `str`
        Full sequence of the protospacer/on-target: 5'-20nt-PAM-3',
        possibly preceded by the upstream sequence.
    offtarget_seq : `str`
        Full sequence of the (off-)target: 5'-20nt-PAM-3', possibly
        preceded by the upstream sequence. If provided, overrules
        the ``mutations`` parameter.

    Returns
    -------
    dna_opening_energy : `tuple` [`float`]
        Free energies required to open the DNA duplex
    rna_duplex_energy : `tuple` [`float`]
        Free energies required to form the RNA duplex
    """


    # Prepare target DNA and guide RNA
    hybrid = GuideTargetHybrid.from_cas9_offtarget(offtarget_seq,
                                                   protospacer)

    # prepare NearestNeighborModel
    nnmodel = NearestNeighborModel
    nnmodel.load_data()
    nnmodel.set_energy_unit("kbt")

    # do calculations
    dna_opening_energy = nnmodel.dna_opening_energy(hybrid)
    rna_duplex_energy = nnmodel.rna_duplex_energy(hybrid)

    # hashable output
    return tuple(dna_opening_energy), tuple(rna_duplex_energy)


def find_average_mm_penalties(protospacer: str,
                              weight: Union[float,
                                            Tuple[float, float]] = None):
    """Find the effective penalties for single point mutations.

    Finds the effective penalties for all possible single point mutations
    on a target, and averages over them to return the position-dependent
    mismatch penalty due to undetermined mismatches.
    """

    # prepare NearestNeighborModel
    nnmodel = NearestNeighborModel
    nnmodel.load_data()
    nnmodel.set_energy_unit("kbt")

    on_target_hybrid = GuideTargetHybrid.from_cas9_protospacer(protospacer)
    u_ontarget = nnmodel.get_hybridization_energy(on_target_hybrid,
                                                  weight=weight)[-1]

    basepairs = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    avg_smm_penalties = np.zeros(20)
    for i in range(20):
        mm_hybr_energies = []

        for nt in basepairs.keys():
            ps_nt = on_target_hybrid.target.seq2[-(i+1)]
            if nt == ps_nt:
                continue

            off_target_hybrid = on_target_hybrid.apply_point_mut(
                f"{ps_nt}{i+1:02d}{nt}"  # e.g. A04T
            )
            u_final = nnmodel.get_hybridization_energy(off_target_hybrid,
                                                       weight=weight)[-1]
            mm_hybr_energies += [u_final]

        avg_smm_penalties[i] = np.mean(mm_hybr_energies) - u_ontarget

    return avg_smm_penalties


class TargetDna:
    """
    Double-stranded DNA site to be opened during R-loop formation.

    Attributes
    ----------
    seq2 : `str`
        The "target sequence", as present on the nontarget DNA strand
        (=protospacer), in 5'-to-3'notation.
    seq1 : `str`
        The target strand (=spacer), in 3'-to-5' notation
    upstream_bp: `str`
        The basepair upstream (5'-side) of the nontarget strand.
    dnstream_bp: `str`
        The basepair downstream (3'-side) of the nontarget strand. For Cas9,
        corresponds to the last basepair of the PAM.
    """

    bp_map = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

    def __init__(self, target_sequence,
                 upstream_nt: str = None,
                 downstream_nt: str = None):

        # non-target strand
        self.seq2 = target_sequence
        # target strand
        self.seq1 = self.__reverse_transcript(self.seq2)

        self.upstream_nt = upstream_nt
        self.downstream_nt = downstream_nt

        if upstream_nt is None:
            self.upstream_bp = None
        else:
            self.upstream_bp = (upstream_nt + "-" +
                                self.bp_map[upstream_nt])
        if downstream_nt is None:
            self.dnstream_bp = None
        else:
            self.dnstream_bp = (downstream_nt + "-" +
                                self.bp_map[downstream_nt])

    def __str__(self):
        """Generate a handy string representation of the DNA duplex."""

        strand2 = self.seq2
        strand1 = self.seq1

        if self.upstream_bp:
            strand2 = self.upstream_bp[0] + strand2
            strand1 = self.upstream_bp[-1] + strand1
        else:
            strand2 = "N" + strand2
            strand1 = "N" + strand1

        if self.dnstream_bp:
            strand2 = strand2 + self.dnstream_bp[0]
            strand1 = strand1 + self.dnstream_bp[-1]
        else:
            strand2 = strand2 + "N"
            strand1 = strand1 + "N"

        # Formatting labels to indicate bp index
        count = 5
        reps = (len(self.seq1) - 3) // count
        labs = ''.join(["{0:<{1}d}".format(count * i, count)
                        for i in reversed(range(1, reps + 1))])[:-1] + '1'
        labs = labs.rjust(len(self.seq1) + 4)
        labs = 4 * ' ' + str(len(self.seq1)).ljust(2) + labs[6:]

        return "\n".join([
            f"3\'-{strand1}-5\' (DNA TS)",
            3 * " " + (2 + len(self.seq1)) * "|",
            f"5\'-{strand2}-3\' (DNA NTS)",
            labs
        ])

    @classmethod
    def __reverse_transcript(cls, sequence: str) -> str:
        """Gives complementary sequence (in opposite direction!)"""
        return ''.join([cls.bp_map[n] for n in sequence])

    @classmethod
    def from_cas9_target(cls, full_target: str) -> 'TargetDna':
        """Make a TargetDna instance from a cas9 target sequence string.

        Assumes a length of 20 nt, followed by a PAM. If the target
        sequence is longer than 23, it includes the upstream nt as well.
        """

        # no upstream nt, target + PAM
        if full_target[-2:] != "GG":
            raise ValueError("Full target should end with 5'-NGG-3' PAM.")

        downstream_nt = full_target[-3]
        upstream_nt = (None if len(full_target) <= 23
                       else full_target[-24])
        target_seq = full_target[-23:-3]
        return cls(target_seq, upstream_nt, downstream_nt)

    def apply_point_mut(self, mutation: str):
        """Change DNA hybrid according to a single point mutation.

        Mutation strings have the form A02T, where the NTS nucleotide A
        at position 2 would get replaced by a nucleotide T.
        """

        old_nt = mutation[0]
        pmut_pos = len(self.seq2) - int(mutation[1:3])
        new_nt = mutation[3]

        if self.seq2[pmut_pos] != old_nt:
            raise ValueError("Mutation doesn't match target sequence")

        new_seq = self.seq2[:pmut_pos] + new_nt
        if pmut_pos < len(self.seq2):
            new_seq += self.seq2[pmut_pos+1:]

        return type(self)(
            new_seq, upstream_nt=self.upstream_nt,
            downstream_nt=self.downstream_nt,
        )

    @classmethod
    def make_random(cls, length: int, seed=None) -> 'TargetDna':
        """Make a random target dna of specified length."""
        random.seed(seed)
        nucleotides = list(cls.bp_map.keys())
        seq = ''.join(random.choices(nucleotides, k=length + 2))
        return cls(upstream_nt=seq[0],
                   target_sequence=seq[1:length+1],
                   downstream_nt=seq[-1])


class GuideTargetHybrid:
    """A ssRNA guide interacting with ds DNA site through R-loop formation.

    Attributes
    ----------
    guide : `str`
        The RNA guide strand, in 5'-to-3' notation
    target : `TargetDna`
        The dsDNA site to be interrogated
    state : `int`
        Length of the R-loop. Only for illustration purposes for now.
    """

    bp_map = {'A': 'U', 'C': 'G', 'G': 'C', 'T': 'A'}

    def __init__(self, guide: str, target: TargetDna, state: int = 0):
        self.guide = guide
        self.target = target
        self.state = state

    def __str__(self):
        """Generate a handy string representation of the R-loop."""
        dna_repr = str(self.target).split('\n')
        hybrid_bps = (''.join(map(str, self.find_mismatches()))
                      .replace('1', "\u00B7")
                      .replace('0', '|'))[::-1]
        hybrid_bps = (hybrid_bps[-self.state:]).rjust(4 + len(self.guide))
        if self.state == 0:
            hybrid_bps = ''

        return "\n".join([
            f" 5\'-{self.guide}-3\'  (RNA guide)",
            hybrid_bps,
            dna_repr[0],
            dna_repr[1][:4 + len(self.guide) - self.state] +
            self.state * " " + "|",
            dna_repr[2],
            dna_repr[3]
        ])

    @classmethod
    def from_cas9_protospacer(cls, protospacer: str, mismatches: str = '',
                              state: int = 0) -> 'GuideTargetHybrid':
        """Instantiate from protospacer and point mutations.

        Parameters
        ----------
        protospacer : `str`
            Full sequence of the protospacer/on-target: 5'-20nt-PAM-3',
            possibly preceded by the upstream sequence.
        mismatches : `str`
            Mismatch desciptors (in the form "A02T") describing how the
            target deviates from the ``protospacer``. Multiple mismatches
            should be space-separated.
        state : `int`
            R-loop hybridization state
        """

        ontarget = TargetDna.from_cas9_target(protospacer)
        guide_rna = ''.join([cls.bp_map[n] for n in ontarget.seq1])

        target = ontarget
        if mismatches != '':
            for mm in mismatches.split():
                target = target.apply_point_mut(mm)

        return cls(guide_rna, target, state)

    @classmethod
    def from_cas9_offtarget(cls, offtarget_seq: str, protospacer: str,
                            state: int = 0) -> 'GuideTargetHybrid':
        """Instantiate from protospacer and point mutations.

        Parameters
        ----------
        offtarget_seq : `str`
            Full sequence of the (off-)target: 5'-20nt-PAM-3', possibly
            preceded by the upstream sequence.
        protospacer : `str`
            Full sequence of the protospacer/on-target: 5'-20nt-PAM-3',
            possibly preceded by the upstream sequence.
        state : `int`
            R-loop hybridization state
        """

        ontarget = TargetDna.from_cas9_target(protospacer)
        guide_rna = ''.join([cls.bp_map[n] for n in ontarget.seq1])
        offtarget = TargetDna.from_cas9_target(offtarget_seq)
        return cls(guide_rna, offtarget, state)

    def apply_point_mut(self, mutation: str) -> 'GuideTargetHybrid':
        return type(self)(
            guide=self.guide,
            target=self.target.apply_point_mut(mutation),
            state=self.state
        )

    def set_rloop_state(self, rloop_state):
        self.state = rloop_state

    def find_mismatches(self):
        """Identify the positions of mismatching guide-target basepairs."""
        return list(find_mismatches_cached(self.target.seq1, self.guide))

    def get_mismatch_pattern(self) -> MismatchPattern:
        return MismatchPattern(self.find_mismatches())


@lru_cache
def find_mismatches_cached(seq1, guide):
    """"Identify the positions of mismatching guide-target basepairs (cached)."""
    bp_map = GuideTargetHybrid.bp_map
    mismatches = []
    for i, n in enumerate(reversed(seq1)):
        mismatches += [0 if bp_map[n] == guide[-1 - i]
                       else 1]
    return tuple(mismatches)


class NearestNeighborModel:
    """A model to estimate nucleic acid stability.

    An implementation of the nearest neighbor model predicting energies
    for guide RNA-target DNA R-loops. Instantiating this class is only
    necessary to load the parameter files, a single object can be used
    to make all energy landscapes.

    Attributes
    ----------
    energy_unit : {'kbt', 'kcalmol'}
        Unit of ouput free energy. For kBT, assuming a temperature of 20°C.

    Notes
    -----
    Method adapted from Alkan et al. (2018). DNA duplex parameters from
    SantaLucia & Hicks (2004), RNA-DNA hybrid duplex parameters from
    Alkan et al. (2018).

    There are 4 contributions to the R-loop energy.

    1) Basestacks in the DNA duplex that should be broken. These
        parameters can be loaded directly from the SantaLucia & Hicks
        dataset. Unlike Alkan et al., we also consider basestacks with the
        basepairs flanking the target region. If these are unknown,
        we take the average energy from all 4 possible basestacks.
    2) Basestacks in the RNA/DNA hybrid that are created. Some of these
        energies are experimentally determined, others are an average
        of dsDNA and dsRNA values.
    3) Internal loops, corresponding to (regions of) mismatches flanked
        by matching basepairs. For internal loops of length 1 and 2,
        these have specific energies, for length > 2, their energies
        are the sum of the left and right basestack and a
        length-specific energy contribution.
    4) Basepair terminals at the end and beginning of the R-loop.
        Alkan et al. consider only external loops, which appear only when
        the guide-target hybrid starts or ends with a mismatch, but
        we always consider the energy contribution due to the first and
        last matching basepair. These energies are typically quite
        small.

    References
    ----------
    .. [1] Alkan F, Wenzel A, Anthon C, Havgaard JH, Gorodkin J (2018).
        CRISPR-Cas9 off-targeting assessment with nucleic acid duplex
        energy parameters. doi.org/10.1186/s13059-018-1534-x
    .. [2] SantaLucia J, Hicks D (2004). The Thermodynamics of DNA
        Structural Motifs. doi.org/10.1146/annurev.biophys.32.110601.141800
    """

    # paths relative to crisprzip source root
    dna_dna_params_file = "data/nucleicacid/santaluciahicks2004.json"
    rna_dna_params_file = "data/nucleicacid/alkan2018.json"
    dna_dna_params: dict = None
    rna_dna_params: dict = None

    energy_unit = "kbt"  # alternative: kcalmol

    @classmethod
    def load_data(cls, force=False):

        # find git root
        def find_git_root(path):
            while not os.path.isdir(os.path.join(path, '.git')):
                parent = os.path.dirname(path)
                if parent == path:
                    raise FileNotFoundError("No .git directory found")
                path = parent
            return path

        # navigate to git root
        current_dir = Path().resolve()
        package_root = find_git_root(current_dir)

        if cls.dna_dna_params is None or force:
            with open(package_root.joinpath(cls.dna_dna_params_file), 'rb') as file:
                cls.dna_dna_params = json.load(file)

        if cls.rna_dna_params is None or force:
            with open(package_root.joinpath(cls.rna_dna_params_file), 'rb') as file:
                cls.rna_dna_params = json.load(file)
    
        with open(package_root.joinpath(cls.dna_dna_params_file), 'rb') as file:
            cls.dna_dna_params = json.load(file)

        with open(package_root.joinpath(cls.rna_dna_params_file), 'rb') as file:
            cls.rna_dna_params = json.load(file)

    @classmethod
    def set_energy_unit(cls, unit: str):
        cls.energy_unit = unit

    @classmethod
    def convert_units(cls, energy_value: Union[float, np.ndarray]):

        if cls.energy_unit == "kcalmol":
            return energy_value

        elif cls.energy_unit == "kbt":
            ref_temp = 310.15  # 310.15 deg K = 37 deg C
            gas_constant = 1.9872E-3  # R = N_A * k_B [units kcal / (K mol)]
            return energy_value / (gas_constant * ref_temp)

    @classmethod
    def get_hybridization_energy(cls, hybrid: GuideTargetHybrid,
                                 weight: Union[float,
                                               Tuple[float, float]] = None) \
            -> np.ndarray:
        """Calculate the R-loop cost.

        Calculates theenergy that is required to open an R-loop
        between the guide RNA and target DNA of the hybrid object
        for each R-loop length. Converts energy units if necessary.

        Parameters
        ----------
        hybrid : `GuideTargetHybrid`
            Hybrid object of which the hybridization energies are calculated
        weight : `float` or `tuple`[`flaot`], optional
            Optional weighing of the dna opening energy and rna duplex energy.
            If `None` (default), no weighing is applied. If `float`, both DNA and
            RNA energies are multiplied by the weight parameter. If `tuple``
            of two `float`s, the first value is used as a multiplier for the
            DNA opening energy, and the second is used as a multiplier for the
            RNA-DNA hybridization energy.

        Returns
        -------
        energy : `numpy.ndarray`
            The energy required for hybridization (in the desired units of
            energy), for each step in the R-loop formation process.
        """

        dna_opening_energy = cls.dna_opening_energy(hybrid)
        rna_duplex_energy = cls.rna_duplex_energy(hybrid)
        if weight is None:
            return cls.convert_units(
                dna_opening_energy + rna_duplex_energy
            )
        elif type(weight) is float:
            return weight * cls.convert_units(
                dna_opening_energy + rna_duplex_energy
            )
        elif type(weight) is tuple:
            return cls.convert_units(
                weight[0] * dna_opening_energy +
                weight[1] * rna_duplex_energy
            )
        else:
            raise ValueError(f"Cannot interpret weight of type {type(weight)}")

    @classmethod
    def dna_opening_energy(cls, hybrid: GuideTargetHybrid) -> np.ndarray:
        """Get the energy required to open the DNA duplex.

        Calculated following the methods from Alkan et al. (2018).
        The DNA opening energy is the sum of all the basestack energies
        in the sequence (negative).
        """

        stacking_energies = cls.dna_dna_params["stacking energies"]

        def average_basestack(side='downstream'):
            # If down- or upstream basepair is unknown, this function
            # loops over all options to find the average basestack
            # energy.
            basestack_energy = 0
            possible_bps = [('A', 'T'), ('C', 'G'),
                            ('G', 'C'), ('T', 'A')]

            for bp in possible_bps:
                if side == 'upstream':
                    bstack = (
                        f"d{bp[0]}{hybrid.target.seq2[0]}/"
                        f"d{bp[1]}{hybrid.target.seq1[0]}"
                    )
                elif side == 'downstream':
                    bstack = (
                        f"d{hybrid.target.seq2[-1]}{bp[0]}/"
                        f"d{hybrid.target.seq1[-1]}{bp[1]}"
                    )
                else:
                    raise ValueError("Side can only be downstream or"
                                     "upstream.")
                basestack_energy += (stacking_energies[bstack] / 4)
            return basestack_energy

        open_energy = np.zeros(len(hybrid.guide) + 1)

        # Rightmost basestack
        if hybrid.target.dnstream_bp is not None:
            basestack = (f"d{hybrid.target.seq2[-1]}"
                         f"{hybrid.target.dnstream_bp[0]}/"
                         f"d{hybrid.target.seq1[-1]}"
                         f"{hybrid.target.dnstream_bp[-1]}")
            open_energy[1:] -= stacking_energies[basestack]
        else:
            open_energy[1:] -= average_basestack("downstream")

        # Handling middle basestacks
        rev_ntstr = hybrid.target.seq2[::-1]
        rev_tgstr = hybrid.target.seq1[::-1]
        for i in (range(0, len(hybrid.guide) - 1)):
            # sorry for the hacky way of formatting basestacks
            basestack = (f"d{rev_ntstr[i:i + 2][::-1]}/"
                         f"d{rev_tgstr[i:i + 2][::-1]}")
            open_energy[i + 1:] -= stacking_energies[basestack]

        # Leftmost basestack
        if hybrid.target.upstream_bp is not None:
            basestack = (f"d{hybrid.target.upstream_bp[0]}"
                         f"{hybrid.target.seq2[0]}/"
                         f"d{hybrid.target.upstream_bp[-1]}"
                         f"{hybrid.target.seq1[0]}")
            open_energy[-1] -= stacking_energies[basestack]
        else:
            open_energy[-1] -= average_basestack("upstream")

        return open_energy

    @classmethod
    def rna_duplex_energy(cls, hybrid: GuideTargetHybrid) -> np.ndarray:
        """Get the energy required to create the RNA:DNA duplex.

        Calculated following the methods from Alkan et al. (2018).
        The RNA duplex energy has three contributions: 1) basestacks,
        2) internal loops, 3) external loops / terminals. Alkan et al.
        only look at external loops, but here, we instead look
        at both basepair terminals, whether or not they are part of
        an external loop."""

        loop_energies = cls.rna_dna_params['loop energies']
        terminal_energies = cls.rna_dna_params['terminal penalties']
        stacking_energies = cls.rna_dna_params['stacking energies']

        def basestack_energy(i: int):
            # Locates basestacks in R-loop of length i and
            # returns their total energy.
            if i < 2:
                return 0.

            energy = 0.
            mm_pos = hybrid.find_mismatches()[:i]
            partial_rev_guide = hybrid.guide[-1:-1 - i:-1]
            partial_rev_tgstr = hybrid.target.seq1[-1:-1 - i:-1]

            for j, nt in enumerate(partial_rev_guide[:-1]):
                if mm_pos[j] == 0 and mm_pos[j + 1] == 0:
                    basestack = f"r{partial_rev_guide[j + 1]}" \
                                f"{partial_rev_guide[j]}/" \
                                f"d{partial_rev_tgstr[j + 1]}" \
                                f"{partial_rev_tgstr[j]}"
                    energy += stacking_energies["2mer"][basestack]
            return energy

        def internal_loops_energy(i: int):
            # Locates internal loops in R-loop of length i and
            # returns their total energy.
            if i < 3:
                return 0.

            energy = 0.
            mm_pos = hybrid.find_mismatches()[:i]
            partial_rev_guide = hybrid.guide[-1:-1 - i:-1]
            partial_rev_tgstr = hybrid.target.seq1[-1:-1 - i:-1]

            # scanning hybrid from right to left
            for j in range(1, i - 1):

                # identify right side of internal loop
                if mm_pos[j - 1] == 0 and mm_pos[j] == 1:
                    looplen = 1

                    # find left side of internal loop
                    for k in range(j + 1, i):
                        if mm_pos[k] == 0:
                            break
                        looplen += 1
                    # if no closing basepair has been found; stop
                    if j + looplen == i:
                        break

                    if looplen == 1:
                        basestacks = (
                            f"r{partial_rev_guide[j + 1]}"
                            f"{partial_rev_guide[j]}"
                            f"{partial_rev_guide[j - 1]}/"
                            f"d{partial_rev_tgstr[j + 1]}"
                            f"{partial_rev_tgstr[j]}"
                            f"{partial_rev_tgstr[j - 1]}"
                        )
                        energy += stacking_energies["3mer"][basestacks]

                    elif looplen == 2:
                        basestacks = (
                            f"r{partial_rev_guide[j + 2]}"
                            f"{partial_rev_guide[j + 1]}"
                            f"{partial_rev_guide[j]}"
                            f"{partial_rev_guide[j - 1]}/"
                            f"d{partial_rev_tgstr[j + 2]}"
                            f"{partial_rev_tgstr[j + 1]}"
                            f"{partial_rev_tgstr[j]}"
                            f"{partial_rev_tgstr[j - 1]}"
                        )
                        energy += stacking_energies["4mer"][basestacks]

                    elif looplen >= 3:
                        left_basestack = (
                            f"r{partial_rev_guide[j + looplen]}"
                            f"{partial_rev_guide[j + looplen - 1]}/"
                            f"d{partial_rev_tgstr[j + looplen]}"
                            f"{partial_rev_tgstr[j + looplen - 1]}"
                        )
                        right_basestack = (
                            f"r{partial_rev_guide[j]}"
                            f"{partial_rev_guide[j - 1]}/"
                            f"d{partial_rev_tgstr[j]}"
                            f"{partial_rev_tgstr[j - 1]}"
                        )
                        if looplen <= 9:
                            loop_energy = loop_energies[f"{2 * looplen} nt"]
                        else:
                            loop_energy = 4.5
                        energy += (
                                stacking_energies["2mer"][left_basestack] +
                                stacking_energies["2mer"][right_basestack] +
                                loop_energy
                        )

            return energy

        def external_loops_energy(i: int):
            # Locates external loops in R-loop of length i and
            # returns their total energy. Unused right now because
            # we're considering all terminal energies instead.

            mm_pos = hybrid.find_mismatches()[:i]
            partial_rev_guide = hybrid.guide[-1:-1 - i:-1]
            partial_rev_tgstr = hybrid.target.seq1[-1:-1 - i:-1]

            if sum(mm_pos) == len(mm_pos):
                return 0.

            energy = 0.

            # open right side
            if mm_pos[0] == 1:
                # find right side of internal loop
                for k in range(0, i):
                    if mm_pos[k] == 0:
                        basepair = f"r{partial_rev_guide[k]}-" \
                                   f"d{partial_rev_tgstr[k]}"
                        energy += terminal_energies[basepair]
                        break

            # open left side
            if mm_pos[-1] == 1:
                # find right side of internal loop
                for k in range(i - 1, 0, -1):
                    if mm_pos[k] == 0:
                        basepair = f"r{partial_rev_guide[k]}-" \
                                   f"d{partial_rev_tgstr[k]}"
                        energy += terminal_energies[basepair]
                        break

            return energy

        def terminals_energy(i: int):
            # Locates terminals in R-loop of length i and
            # returns their total energy.

            mm_pos = hybrid.find_mismatches()[:i]
            partial_rev_guide = hybrid.guide[-1:-1 - i:-1]
            partial_rev_tgstr = hybrid.target.seq1[-1:-1 - i:-1]

            if sum(mm_pos) == len(mm_pos):
                return 0.

            energy = 0.
            # find rightmost basepair
            for k in range(0, i):
                if mm_pos[k] == 0:
                    basepair = f"r{partial_rev_guide[k]}-" \
                               f"d{partial_rev_tgstr[k]}"
                    energy += terminal_energies[basepair]
                    break

            # find leftmost basepair
            for k in list(np.arange(i, 0, -1) - 1):
                if mm_pos[k] == 0:
                    basepair = f"r{partial_rev_guide[k]}-" \
                               f"d{partial_rev_tgstr[k]}"
                    energy += terminal_energies[basepair]
                    break

            return energy

        # Looping over R-loop states, adding all energy contributions
        duplex_energy = np.array([(
            basestack_energy(i) +
            internal_loops_energy(i) +
            terminals_energy(i)
        ) for i in range(len(hybrid.guide) + 1)
        ])
        return duplex_energy
