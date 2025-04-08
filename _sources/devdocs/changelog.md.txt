# Changelog

All notable changes to CRISPRzip will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1.1] - 2025-04-08
### Bugfix
- More robust approach to checking sequence input. Any sequence (protospacer/off-target)
  can now be offered in a 20/23/24 nt format.

## [1.1.0] - 2025-02-18
### Added
- Function `load_landscape` in `crisprzip.kinetics` for robust parameter loading
- Cross-platform testing: MacOS, Ubuntu, Windows for Python 3.10, 3.11, 3.12

### Bugfix
- Parameters for nucleic acid stability and landscape definitions are now included as package data which avoids relative path issues.


## [1.0.0] - 2025-02-04
### Added
- Initial release of the package.
- Comprehensive documentation for setup and usage, hosted on GitHub-pages.
- Implemented end-to-end testing.

### Notes
- This version marks the first stable release of CRISPRzip.
- Feedback and contributions are welcome to improve the project.
