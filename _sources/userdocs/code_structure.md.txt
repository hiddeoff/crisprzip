# Code structure

## Modules
The code for CRISPRzip is set up according to OOP principles. Often, you first 
have to instantiate an object (like `Searcher` or `SearcherPlotter`), after 
which you can use its methods to do calculations etcetera.

In this setup, there is a strict hierarchy: if class B has any reference to 
class A, then class A should never reference class B. This design principle 
allows easy redefinition of higher-level classes like class B. Another feature 
of OOP that we have made us of is inheritance. For instance, all functionality 
around DNA sequences is built on top of the sequence-agnostic model, by having
child classes with just one or two methods different than the parent.

(sec-kinetics)=
### crisprzip.kinetics
🔗 *go to [API reference](../apidocs/crisprzip.kinetics.rst)*

The `kinetics` module contains the core functionality of CRISPRzip. Its classes
define how RNP landscapes should be generated. There are classes that are
independent of sequence,
- `Searcher`;
- `SearcherTargetComplex`,

and classes that do depend on target sequence,
- `BareSearcher`;
- `GuidedSearcher`;
- `SearcherSequenceComplex`.

To start with the simpler sequence-independent model, the `Searcher` object
corresponds to a RNP (such as Cas9) that has not bound a target yet. It contains
the necessary parameters to construct R-loop landscapes for any target. 
For a particular binding candidate and its according `MismatchPattern`, 
one can make a `SearcherTargetComplex`. With the added information where the 
mismatches are located in the R-loop, the landscape for this target can be 
generated and kinetics can be obtained. For details on 
the landscape generation and how hybridization kinetics follow from it, see 
[Eslami-Mossalam et al. (2022)](#ref-eslami2022).
To make predictions, the Master Equation for this system needs to be solved. 
The implementation of this linear algebra problem can be found in
[crisprzip.matrix_expon](#sec-matrixexpon).

The sequence-dependent classes use the code in [crisprzip.nucleic_acid](#sec-nucleicacid) to
find the cost of a bare R-loop. Hence, their attributes describing the on-target
landscape and mismatch penalties no longer include the cost due to the 
nucleic acid stability, but only represent the sequence-independent contributions 
(e.g. internal stability of the protein and interactions with the DNA backbone).
With slightly different methods that include the cost of an R-loop, it makes landscapes 
that depend on sequence, and accordingly, predicts kinetics that vary with sequence.

(sec-matrixexpon)=
### crisprzip.matrix_expon
🔗 *go to [API reference](../apidocs/crisprzip.matrix_expon.rst)*

To obtain kinetics from a particular landscape, one needs to solve the
Master Equation, which involves exponentiating a 22x22 matrix. As this is a
somewhat costly operation (which is repeated many times when training the
model), we have optimized this code by using the
[Numba](https://numba.readthedocs.io/en/stable/) package. Numba compiles Python 
code 'just-in-time' (upon the first call), resulting in much faster execution.

In this module, there are few versions of the same function. They differ in 1) 
which arguments they can sweep over and 2) whether or not they use Numba.

(sec-nucleicacid)=
### crisprzip.nucleic_acid
🔗 *go to [API reference](../apidocs/crisprzip.nucleic_acid.rst)*

To be able to make sequence-dependent predictions, we need a model for the energetic 
cost of an R-loop. The Nearest-Neighbor model provides us with such a model:
it approximates the stability of a nucleic acid construct by the sum over
all its basestacks (=neighboring basepairs). We adopt an implementation by 
[Alkan et al. (2018)](#ref-alkan2018) for  matching and mismatch stretches of 
DNA:DNA and RNA:DNA. 

(sec-coarsegrain)=
### crisprzip.coarsegrain
🔗 *go to [API reference](../apidocs/crisprzip.coarsegrain.rst)*

By coarse-graining the landscapes that are generated in 
[crisprzip.kinetics](#sec-kinetics), the main features of the landscapes
are shown more clearly. This helps when comparing different landscapes with each 
other or when comparing CRISPRzip rates to smFRET data. Again, the coarse-
graining procedure is explained in [Eslami-Mossalam et al. (2022)](#ref-eslami2022).

(sec-plotting)=
### crisprzip.plotting
🔗 *go to [API reference](../apidocs/crisprzip.plotting.rst)*

With the `SearcherPlotter` object, you can visualize the attributes of 
a `Searcher` object, like its on- or off-target landscapes.


### References
(ref-eslami2022)=
Eslami-Mossallam B et al. (2022) *A kinetic model predicts SpCas9 activity,
improves off-target classification, and reveals the physical basis of
targeting fidelity.* Nature Communications.
[10.1038/s41467-022-28994-2](https://doi.org/10.1038/s41467-022-28994-2)

(ref-alkan2018)=
Alkan et al. (2018) *CRISPR-Cas9 off-targeting assessment with nucleic acid 
duplex energy parameters.* Genome Biology.
[10.1186/s13059-018-1534-x](https://doi.org/10.1186/s13059-018-1534-x)
