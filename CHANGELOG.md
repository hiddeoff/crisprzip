# Changelog

All notable changes to CRISPRzip will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.2.1] - 2025-07-31
### Fixed
- Sequence and average parameter sets are updated to correspond with
  the proper R-loop cost calculations.

### Added
- Functionality to save R-loop cost weighing as parameter 'weight' to
  parameter sets and load them as a BareSearcher attribute.
- Function to load binding rates from parameter set files.


## [1.2.0] - 2025-07-08
### Fixed
- Major bugfix: R-loop cost is now properly converted to units kBT (at 
  room temperature, 20 Celsius). In all previous versions, R-loop cost 
  was reported in units kcal/mol, being a factor ~1.7 too small.
- End-to-end test no longer passes automatically. Test values are updated 
  according to new R-loop costs.

### Added
- Inclusion of Nearest-Neighbor values from Banerjee et al. (2020) for 
  RNA/DNA hybrids at physiological salt conditions (100 mM NaCl). Support 
  for these values is quite limited; they are only available for matching
  hybrid basepairs, and no similar salt correction can yet be made for 
  the dsDNA base pairs.
- A set of sequence parameters for a smooth protein landscape. However, 
  these values will be replaced soon (see Notes below).

### Notes
- With the R-loop cost being corrected, new sequence-dependent energy 
  landscapes should be generated and included with this package. This makes
  version 1.2.0 a temporary version, and the more stable 1.2.1 will
  follow very soon.


## [1.1.1] - 2025-04-08
### Fixed
- More robust approach to checking sequence input. Any sequence (protospacer/off-target)
  can now be offered in a 20/23/24 nt format.

## [1.1.0] - 2025-02-18
### Added
- Function `load_landscape` in `crisprzip.kinetics` for robust parameter loading
- Cross-platform testing: MacOS, Ubuntu, Windows for Python 3.10, 3.11, 3.12

### Fixed
- Parameters for nucleic acid stability and landscape definitions are now included as package data which avoids relative path issues.


## [1.0.0] - 2025-02-04
### Added
- Initial release of the package.
- Comprehensive documentation for setup and usage, hosted on GitHub-pages.
- Implemented end-to-end testing.

### Notes
- This version marks the first stable release of CRISPRzip.
- Feedback and contributions are welcome to improve the project.
