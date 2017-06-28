# Change Log
All notable changes to this project will be documented in this file.

## [v0.1.2]
### Fixed
- Fixed a typo in the MANIFEST.in file that caused the header files
  not to be distributed on PyPI.

## [v0.1.1]
### Added
- Additional fields for PyPI support were added to the setup.py file.

### Fixed
- Fixed typos in the version numbers in various files.

## [v0.1.0]
### Added
- A method `getCppWLC2D` was added to helpers.py for generating
  two-dimensional wormlike chains.
		
## [v0.0.1]
### Fixed
- The Python interface to the self-avoiding wormlike chain no longer
  throws an error when the input path length is less than three.

## v0.0.0
### Added
- Documentation
- Some minimal unit tests and doctests
- Verification routine to check the accuracy of the WLC algorithms

### Changed
- The project was restructured to better match a standard Python/C++
  open-source format.

### Fixed
- Fixed off-by-one error in the wormlike chain generation code.
- Errors in computation of chain statistics

[v0.1.2]: https://github.com/kmdouglass/PolymerCpp/compare/v0.1.1...v0.1.2
[v0.1.1]: https://github.com/kmdouglass/PolymerCpp/compare/v0.1.0...v0.1.1
[v0.1.0]: https://github.com/kmdouglass/PolymerCpp/compare/v0.0.1...v0.1.0
[v0.0.1]: https://github.com/kmdouglass/PolymerCpp/compare/v0.0.0...v0.0.1
