"""Simulate the kinetics of R-loop formation by a CRISPR(-like)
endonuclease."""

import importlib.resources

import numpy as np

from .matrix_expon import *
from .nucleic_acid import *


class Searcher:
    """A (CRISPR-associated) RNA-guided endonuclease.

    Attributes
    ----------
    guide_length : `int`
        Length of the nucleic acid guide N (in bp)
    on_target_landscape : `numpy.ndarray`, (N,)
        Contains the hybridization energies of the intermediate R-loop
        states on an on-target, relative to the PAM energy.
    mismatch_penalties : `numpy.ndarray`, (N,)
        Contains the energetic penalties associated with a mismatch
        at a particular R-loop position.
    internal_rates : `dict` [`str`, `float`]
        Specifies the internal (=context-independent) rates in the
        model:

        ``"k_off"``
            PAM-unbinding rate in s⁻¹ (`float`).
        ``"k_f"``
            Forward rate in s⁻¹ (`float`).
        ``"k_clv"``
            Cleavage rate in s⁻¹ (`float`).
    pam_detection : `bool`, optional
        If true, the landscape includes a PAM recognition state. `True`
        by default.
    """

    def __init__(self,
                 on_target_landscape: ArrayLike,
                 mismatch_penalties: ArrayLike,
                 internal_rates: dict,
                 pam_detection=True, *args, **kwargs):

        # convert on_target_landscape and mismatch_penalties to numpy
        if (type(on_target_landscape) != np.ndarray or
                type(mismatch_penalties) != np.ndarray):
            on_target_landscape = np.array(on_target_landscape)
            mismatch_penalties = np.array(mismatch_penalties)

        # check whether parameters are 1d arrays
        if on_target_landscape.ndim > 1 or mismatch_penalties.ndim > 1:
            raise ValueError('Landscape parameters must be 1d arrays')

        # check whether landscape dimensions agree with guide length
        guide_length = mismatch_penalties.size
        if on_target_landscape.size != pam_detection + guide_length - 1:
            raise ValueError('Landscape dimensions do not match guide length')

        # check whether internal_rates dictionary contains proper keys
        if not ('k_off' in internal_rates and
                'k_f' in internal_rates and
                'k_clv' in internal_rates):
            raise ValueError('Forward rates dictionary should include k_off, '
                             'k_f and k_clv as keys')

        # assign values
        self.guide_length = guide_length
        self.pam_detection = pam_detection

        self.on_target_landscape = on_target_landscape
        self.mismatch_penalties = mismatch_penalties
        self.internal_rates = internal_rates

    def get_forward_rate_array(self, k_on: float, dead=False) -> np.ndarray:
        """Turn the forward rate dictionary into proper array."""
        if not dead:
            k_clv = self.internal_rates['k_clv'] * np.ones(1)
        else:
            k_clv = np.zeros(1)

        forward_rate_array = np.concatenate(
            #  solution state
            (k_on * np.ones(1),
             # PAM and intermediate R-loop states
             self.internal_rates['k_f'] *
             np.ones(self.on_target_landscape.size),
             # final/complete R-loop state
             k_clv,
             # cleaved state
             np.zeros(1))
        )
        return forward_rate_array

    def probe_target(self, target_mismatches: MismatchPattern) \
            -> 'SearcherTargetComplex':
        """Form a SearcherTargetComplex by binding a target with
        mismatch pattern ``target_mismatches``."""
        return SearcherTargetComplex(self.on_target_landscape,
                                     self.mismatch_penalties,
                                     self.internal_rates,
                                     target_mismatches)

    def calculate_solution_energy(self, k_on):
        """Given an on-rate, returns the effective free energy of the
        solution state (under the assumption of detailed balance)"""
        return np.log(k_on / self.internal_rates['k_off'])


class BareSearcher(Searcher):
    """`Searcher` with only protein contribution to the energy landscapes.

    The `BareSearcher` has no built-in energy contributions from
    nucleic acid; its energy parameters are only due to the nonspecific
    interactions between the protein and the target DNA. The difference
    between the `Searcher` and `BareSearcher` values of the
    ``on_target_landscape`` and ``mismatch_penalties`` attributes are
    due to nucleic acid interactions as determined in the nucleic_acid module.
    """

    @classmethod
    def from_searcher(cls, searcher, protospacer: str,
                      weight: Union[float, Tuple[float, float]] = None,
                      *args, **kwargs) -> 'BareSearcher':
        """Generate object from default `Searcher`

        Subtracts the nearest-neighbour hybridization energies from
        the definition of a 'normal' `Searcher`. Averages over the
        penalties due to single point mutations to find the difference
        with the original ``mm_penalties`` variable.

        Parameters
        ----------
        searcher : `Searcher`
            Object to be transformed.
        protospacer : `str`
            Full sequence of the protospacer/on-target. Can be provided in 3 formats:

            - 20 nts: 5'-target-3'. All nucleotides should be specified.
            - 23 nts: 5'-target-PAM-3'. The PAM should be specified or
                provided as 'NGG'.
            - 24 nts: 5'-upstream_nt-target-PAM-3'. The upstream_nt can be
                specified or provided as 'N'.

        weight : `float` or `tuple` [`float`], optional
            Optional weighing of the DNA opening energy and RNA duplex energy.
            If None (default), no weighing is applied. If float, both DNA and
            RNA energies are multiplied by the weight parameter. If `tuple`
            of two `float`s, the first value is used as a multiplier for the
            DNA opening energy, and the second is used as a multiplier for the
            RNA-DNA hybridization energy.
        *args
            Additional positional arguments, passed on to the constructor.
        **kwargs
            Additional keyword arguments, passed on to the constructor.

        Returns
        -------
        bare_searcher : `BareSearcher`
            New instance of the `BareSearcher` class.
        """

        ontarget_na_energies = get_hybridization_energy(
            protospacer=protospacer,
            weight=weight
        )
        average_mm_penalties = find_average_mm_penalties(
            protospacer=protospacer,
            weight=weight
        )

        return cls(
            on_target_landscape=(searcher.on_target_landscape -
                                 ontarget_na_energies[1:]),
            mismatch_penalties=(searcher.mismatch_penalties -
                                average_mm_penalties),
            internal_rates=searcher.internal_rates,
            pam_detection=searcher.pam_detection,
            protospacer=protospacer,  # for the GuidedSearcher constructor
            *args, **kwargs
        )

    def to_searcher(self, protospacer: str,
                    weight: Union[float,
                                  Tuple[float, float]] = None) -> Searcher:
        """Turn `BareSearcher` into default `Searcher`.

        Add the nearest-neighbour hybridization energies to
        BareSearcher energy parameters to retrieve a 'normal' Searcher
        object with effective landscape and (average) mismatch penalties.

        Parameters
        ----------
        protospacer : `str`
            Full sequence of the protospacer/on-target. Can be provided in 3 formats:

            - 20 nts: 5'-target-3'. All nucleotides should be specified.
            - 23 nts: 5'-target-PAM-3'. The PAM should be specified or
                provided as 'NGG'.
            - 24 nts: 5'-upstream_nt-target-PAM-3'. The upstream_nt can be
                specified or provided as 'N'.

        weight : `float` or `tuple` [`float`], optional
            Optional weighing of the DNA opening energy and RNA duplex energy.
            If None (default), no weighing is applied. If float, both DNA and
            RNA energies are multiplied by the weight parameter. If `tuple`
            of two `float`s, the first value is used as a multiplier for the
            DNA opening energy, and the second is used as a multiplier for the
            RNA-DNA hybridization energy.

        Returns
        -------
        searcher : `Searcher`
            New instance of the `Searcher` class.
        """

        ontarget_na_energies = get_hybridization_energy(
            protospacer=protospacer,
            weight=weight
        )
        average_mm_penalties = find_average_mm_penalties(
            protospacer=protospacer,
            weight=weight
        )

        return Searcher(
            self.on_target_landscape + ontarget_na_energies[1:],
            self.mismatch_penalties + average_mm_penalties,
            self.internal_rates,
            self.pam_detection
        )

    def bind_guide_rna(self, protospacer: str,
                       weight: Union[float, Tuple[float, float]] = None) \
            -> 'GuidedSearcher':
        """Create a GuidedSearcher object."""
        return GuidedSearcher(
            on_target_landscape=self.on_target_landscape,
            mismatch_penalties=self.mismatch_penalties,
            internal_rates=self.internal_rates,
            pam_detection=self.pam_detection,
            protospacer=protospacer,
            weight=weight,
        )

    def probe_target(self, target_mismatches: MismatchPattern) \
            -> 'SearcherTargetComplex':
        """Prohibited, see `BareSearcher.probe_sequence`."""
        raise ValueError("A Bare/GuidedSearcher object cannot probe a "
                         "target without a defined sequence. Use "
                         "probe_sequence() instead.")

    def probe_sequence(self, protospacer: str, target_seq: str,
                       weight: Union[float, Tuple[float, float]] = None) \
            -> 'SearcherSequenceComplex':
        """Instantiate a `SearcherSequenceComplex`."""
        return SearcherSequenceComplex(self.on_target_landscape,
                                       self.mismatch_penalties,
                                       self.internal_rates,
                                       protospacer=protospacer,
                                       target_seq=target_seq,
                                       weight=weight)


class GuidedSearcher(BareSearcher):
    """`BareSearcher` with predefined protospacer/gRNA."""

    def __init__(self, on_target_landscape: np.ndarray,
                 mismatch_penalties: np.ndarray, internal_rates: dict,
                 protospacer: str,
                 weight: Union[float, Tuple[float, float]] = None,
                 *args, **kwargs):
        super().__init__(on_target_landscape=on_target_landscape,
                         mismatch_penalties=mismatch_penalties,
                         internal_rates=internal_rates,
                         *args, **kwargs)
        self.protospacer = None
        self.guide_rna = None
        self.weight = weight
        self.set_on_target(protospacer)

    def to_searcher(self, *args, **kwargs) -> Searcher:
        """Turn `GuidedSearcher` into default `Searcher`."""
        return super().to_searcher(self.protospacer, self.weight)

    def to_bare_searcher(self) -> BareSearcher:
        """Turn `GuidedSearcher` into default `BareSearcher`."""
        return BareSearcher(
            on_target_landscape=self.on_target_landscape,
            mismatch_penalties=self.mismatch_penalties,
            internal_rates=self.internal_rates
        )

    def set_on_target(self, protospacer: str) -> None:
        """Set the on-target sequence."""
        self.protospacer = protospacer
        self.guide_rna = (GuideTargetHybrid
                          .from_cas9_protospacer(protospacer)
                          .guide)

    def probe_sequence(self, target_seq: str,
                       weight: Union[float, Tuple[float, float]] = None,
                       *args, **kwargs) -> 'SearcherSequenceComplex':
        """Instantiate a `SearcherSequenceComplex`."""
        if weight is None:
            weight = self.weight
        return super().probe_sequence(self.protospacer, target_seq,
                                      weight=weight)


class SearcherTargetComplex(Searcher):
    """An RNA-guided endonuclease bound to a particular (off-)target.

    Attributes
    ----------
    target_mismatches : `numpy.ndarray`, (N,)
        Positions of mismatches in the guide-target hybrid, with entries
        0 (matches) and 1 (mismatches).
    off_target_landscape : `numpy.ndarray`, (N,)
        Hybridization energies of the intermediate R-loop states on the
        current off-target.
    """

    def __init__(self, on_target_landscape: np.ndarray,
                 mismatch_penalties: np.ndarray, internal_rates: dict,
                 target_mismatches: MismatchPattern, *args, **kwargs):
        Searcher.__init__(self, on_target_landscape, mismatch_penalties,
                          internal_rates, *args, **kwargs)

        # check dimensions of mismatch position array
        if target_mismatches.length != self.guide_length:
            raise ValueError('Target array should be of same length as guide')
        else:
            self.target_mismatches = target_mismatches

        self.off_target_landscape = self._get_off_target_landscape()
        self.backward_rate_array = self._get_backward_rate_array()

    def _get_off_target_landscape(self) -> np.ndarray:
        """Add penalties due to mismatches to landscape."""
        landscape_penalties = np.cumsum(
            self.target_mismatches.pattern *
            self.mismatch_penalties
        )
        if self.pam_detection:
            off_target_landscape = (self.on_target_landscape +
                                    landscape_penalties)
        else:
            raise ValueError('No support yet for off target landscape'
                             'without PAM detection')

        return off_target_landscape

    def _get_landscape_diff(self) -> np.ndarray:
        """Return the difference between landscape states."""
        if self.pam_detection:
            hybrid_landscape = np.concatenate((
                np.zeros(1),  # preceding zero representing the PAM state
                self.off_target_landscape
            ))
        else:
            # add potential penalty on the first state (= not PAM)
            hybrid_landscape = np.concatenate((
                np.array([
                    (self.target_mismatches.pattern[0] *
                     self.mismatch_penalties[0])
                ]),
                self.off_target_landscape
            ))
        return np.diff(hybrid_landscape, prepend=np.zeros(1))

    def _get_backward_rate_array(self) -> np.ndarray:
        """Obtains backward rates from detailed balance condition."""
        boltzmann_factors = np.exp(self._get_landscape_diff())
        backward_rate_array = np.concatenate(
            #  solution state
            (np.zeros(1),
             # PAM state
             np.array([self.internal_rates['k_off']]),
             # R-loop states
             self.internal_rates['k_f'] * boltzmann_factors[1:],
             # cleaved state
             np.zeros(1))
        )
        return backward_rate_array

    def get_rate_matrix(self, on_rate: float, dead=False) -> np.ndarray:
        """Set up the rate matrix describing the master equation."""

        forward_rates = self.get_forward_rate_array(k_on=on_rate, dead=dead)
        backward_rates = self.backward_rate_array

        diagonal1 = -(forward_rates + backward_rates)
        diagonal2 = backward_rates[1:]
        diagonal3 = forward_rates[:-1]
        rate_matrix = (np.diag(diagonal1, k=0) +
                       np.diag(diagonal2, k=1) +
                       np.diag(diagonal3, k=-1))

        return rate_matrix

    def solve_master_equation(self, initial_condition: np.ndarray,
                              time: Union[float, np.ndarray],
                              on_rate: Union[float, np.ndarray],
                              dead=False, rebinding=True, mode='fast') -> \
            np.ndarray:
        """Calculate the occupancy of the landscape over time.

        Calculates how the occupancy of the landscape states evolves by
        evaluating the master equation. Absorbing states (solution and
        cleaved state) are explicitly incorporated. Can vary either
        time or on_rate but not both.

        Parameters
        ----------
        initial_condition : `numpy.ndarray`, (N+3,)
            Vector showing the initial occupancy on the hybridization
            landscape. Should sum to 1.
        time : `float` or `numpy.ndarray`, (M,)
            Times at which the master equation is evaluated.
        on_rate : `float` or `numpy.ndarray` (M,)
            Rate (Hz) with which the searcher binds the target from solution.
        dead : `bool`, optional
            If `True`, cleavage rate is set to zero to simulate the
            catalytically inactive dCas9 variant.
        rebinding : `bool`, optional
            If `True`, ``on-rate`` is left intact, if `False`, ``on-rate``
            is set to zero and solution state becomes absorbing.
        mode : {'fast', 'iterative'}, optional
            If 'fast' (default), uses Numba implementation to do fast
            matrix exponentiation. If 'iterative', uses the
            (~20x slower) iterative procedure. Whenever the fast
            implementation fails, falls back to the iterative.

        Returns
        -------
        landscape_occupancy : `numpy.ndarray`, (N+3,) or (N+3, M)
            Occupancy of the landscape states for specified ``time`` and
            ``on_rate``.
        """

        # check dimensions initial condition
        no_states = 3 + self.on_target_landscape.size
        if initial_condition.shape[0] != no_states:
            initial_condition = initial_condition.T
            if initial_condition.shape[0] != no_states:
                raise ValueError('Initial condition should be of same '
                                 'length as hybridization landscape.')

        # if rebinding is prohibited, on-rate should be zero
        if not rebinding:
            on_rate = 0.

        # determines whether to sweep time or k_on (not both)
        vary_time = (False if ((not isinstance(time, np.ndarray)) or
                               time.size == 1)
                     else True)
        vary_k_on = (False if ((not isinstance(on_rate, np.ndarray)) or
                               on_rate.size == 1)
                     else True)
        vary_init = (False if initial_condition.ndim == 1 else True)

        # variable time and k_on: no support (yet)
        if (vary_time + vary_k_on + vary_init) > 1:
            raise ValueError("Cannot iterate over both multiple arguments."
                             "Choose either variable time, k_on, or init.")

        # unique time & k_on (or variable init): handle as variable time
        if not (vary_time or vary_k_on):
            vary_time = True

        # variable time
        if vary_time:
            rate_matrix = self.get_rate_matrix(on_rate, dead=dead)

            # trivial case
            if not isinstance(time, np.ndarray) and np.isclose(time, 0.):
                return initial_condition.T

            # making sure that time is a 1d ndarray
            time = np.atleast_1d(time)

            # where the magic happens; evaluating the master equation
            if mode == 'fast':
                exp_matrix = exponentiate_fast(rate_matrix, time)

                # this is a safe alternative for exponentiate_fast
                if exp_matrix is None:
                    exp_matrix = exponentiate_iterative(rate_matrix, time)

            elif mode == 'iterative':
                exp_matrix = exponentiate_iterative(rate_matrix, time)
            else:
                raise ValueError(f'Cannot recognize mode {mode}')

        # variable k_on
        elif vary_k_on:
            # This reference rate matrix will be updated repeatedly
            ref_rate_matrix = self.get_rate_matrix(0., dead=dead)

            # where the magic happens; evaluating the master equation
            if mode == 'fast':
                exp_matrix = exponentiate_fast_var_onrate(
                    ref_rate_matrix, float(time), on_rate
                )

                # this is a safe alternative for exponentiate_fast
                if exp_matrix is None:
                    exp_matrix = exponentiate_iterative_var_onrate(
                        ref_rate_matrix, time, on_rate
                    )

            elif mode == 'iterative':
                exp_matrix = exponentiate_iterative_var_onrate(
                    ref_rate_matrix, time, on_rate
                )
            else:
                raise ValueError(f'Cannot recognize mode {mode}')

        # This case should never be true
        else:
            raise Exception

        # Shared final maths for variable time & on_rate

        # calculate occupancy: P(t) = exp(Mt) P0
        landscape_occupancy = exp_matrix.dot(initial_condition)

        # avoid negative occupancy (if present, these should be tiny)
        landscape_occupancy = np.maximum(landscape_occupancy,
                                         np.zeros(landscape_occupancy.shape))

        # normalize P(t) to correct for rounding errors
        total_occupancy = np.sum(landscape_occupancy, axis=0)

        # recognize unsafe entries for division (zero/nan/inf)
        unsafe = np.any(np.stack((total_occupancy == 0.,
                                  np.isnan(total_occupancy),
                                  np.isinf(total_occupancy)),
                                 axis=1),
                        axis=1)

        # normalize or assign nan
        landscape_occupancy[:, ~unsafe] = (landscape_occupancy[:, ~unsafe] /
                                           total_occupancy[~unsafe])
        landscape_occupancy[:, unsafe] = (landscape_occupancy[:, unsafe] *
                                          float('nan'))

        return np.squeeze(landscape_occupancy.T)

    def get_cleaved_fraction(self, time: Union[float, np.ndarray],
                             on_rate: float) -> np.ndarray:
        """Get the fraction of cleaved targets after a specified time.

        Parameters
        ----------
        time : `float` or `numpy.ndarray`, (M,)
            Times at which the cleaved fraction is calculated
        on_rate : `float`
            Rate (Hz) with which the searcher binds the target from solution.

        Returns
        -------
        cleaved_fraction : `float` or `numpy.ndarray` (M,)
            Fraction of targets that is expected to be cleaved at ``time``.
        """

        unbound_state = np.concatenate(
            (np.ones(1), np.zeros(self.on_target_landscape.size + 2))
        )
        prob_distr = self.solve_master_equation(unbound_state, time,
                                                on_rate)
        cleaved_fraction = prob_distr.T[-1]
        return cleaved_fraction

    def get_bound_fraction(self, time: Union[float, np.ndarray],
                           on_rate: Union[float, np.ndarray],
                           pam_inclusion: float = 1.) -> np.ndarray:
        """Get the fraction of bound targets after a specified time.

        This calculation assuming that searcher is catalytically
        dead/inactive.

        Parameters
        ----------
        time : `float` or `numpy.ndarray`, (M,)
            Times at which the master equation is evaluated.
        on_rate : `float` or `numpy.ndarray` (M,)
            Rate (Hz) with which the searcher binds the target from solution.
        pam_inclusion : `float`, optional
            Contribution of the PAM state to the bound fraction. When 1.0
            (default), all PAM-bound searchers contribute to the bound
            fraction, when 0.0, nothing does.

        Returns
        -------
        cleaved_fraction: `float` or `numpy.ndarray` (M,)
            Fraction of targets that is expected to be bound after ``time``
            and with binding rates ``on_rate``.
        """

        if not 0. <= pam_inclusion <= 1.:
            raise ValueError(f"PAM inclusion should be between 0. (no PAM "
                             f"contribution) and 1. (full PAM contribution) "
                             f"but is {pam_inclusion: .1f}.")

        unbound_state = np.concatenate(
            (np.ones(1), np.zeros(self.on_target_landscape.size + 2))
        )
        prob_distr = self.solve_master_equation(unbound_state, time,
                                                on_rate, dead=True)
        bound_fraction = (1 - prob_distr.T[0] +
                          - prob_distr.T[1] * (1 - pam_inclusion))
        return bound_fraction


class SearcherSequenceComplex(GuidedSearcher, SearcherTargetComplex):
    """An RNA-guided endonuclease bound to a particular (off-)target sequence.

    The `SearcherSequenceComplex` is like a `SearcherTargetComplex`,
    but with sequence-specific nucleic acid contributions. The most
    important difference is how it calculates its ``off_target_landscape``,
    the sum of all energetic contributions: protein landscape, protein
    penalties, DNA opening energy and RNA-DNA hybridization energy.
    """

    def __init__(self, on_target_landscape: np.ndarray,
                 mismatch_penalties: np.ndarray, internal_rates: dict,
                 protospacer: str, target_seq: str,
                 weight: Union[float, Tuple[float, float]] = None,
                 *args, **kwargs):

        # extra (sequence-related) attributes
        self.protospacer = protospacer
        self.target_seq = target_seq
        self.hybrid = GuideTargetHybrid.from_cas9_offtarget(
            protospacer=protospacer,
            offtarget_seq=self.target_seq
        )
        target_mismatches = self.hybrid.get_mismatch_pattern()
        self.weight = weight

        super().__init__(
            on_target_landscape=on_target_landscape,
            mismatch_penalties=mismatch_penalties,
            internal_rates=internal_rates,
            target_mismatches=target_mismatches,
            protospacer=protospacer,
            weight=weight,
            *args, **kwargs
        )

    def _get_off_target_landscape(self):
        """Add R-loop cost to the protein landscape."""
        internal_na_energy = get_hybridization_energy(
            protospacer=self.protospacer,
            offtarget_seq=self.target_seq,
            weight=self.weight
        )[1:]
        protein_na_energy = (
            SearcherTargetComplex._get_off_target_landscape(self)
        )
        return protein_na_energy + internal_na_energy


def load_landscape(parameter_set: str):
    """Load a parameter set describing the landscape energies for a Searcher

    Parameters
    ----------
    parameter_set : `str`
        Specifies which parameter set to load. This can be one of
        the default parameter sets, distributed along with crisprzip:
        - `sequence_params`: for sequence-specific kinetics
        - `average_params`: for average kinetics
        - `average_params_legacy`: for average kinetics according to Eslami-Mossalam et al. (2021)
        Alternatively, one can provide a path to a JSON-file that describes
        the parameter set. See the Notes for the requirements for the structure
        of such a file.

    Returns
    -------
    searcher_obj : `Searcher` or `BareSearcher`
        An instance of the `Searcher` or one of its subclasses.

    Notes
    -----
    JSON files containing parameter sets for Searcher objects should contain
    at least the following keys:
    - `searcher_class`, corresponding to a class in `crisprzip.kinetics`;
    - `param_values`, with (at least) the arguments for object instantiation,
        - `on_target_landscape`
        - `mismatch_penalties`
        - `internal_rates`.

    """

    available_paramsets = [
        file.stem for file
        in importlib.resources.files("crisprzip.landscape_params").iterdir()
        if (file.is_file() and file.suffix == '.json')
    ]

    if parameter_set in available_paramsets:
        with (importlib.resources.files("crisprzip.landscape_params")
              .joinpath(f"{parameter_set}.json").open("r") as file):
            landscape_params = json.load(file)

    elif Path(parameter_set).is_file() and Path(parameter_set).suffix == '.json':
        with open(Path(parameter_set), 'rb') as file:
            landscape_params = json.load(file)

    else:
        raise ValueError(f"Could not find '{parameter_set}, neither as a"
                         f"predefined parameter set nor as a custom JSON-file.")

    searcher_cls = globals()[landscape_params["searcher_class"]]
    searcher_obj = searcher_cls(**landscape_params['param_values'])
    return searcher_obj
