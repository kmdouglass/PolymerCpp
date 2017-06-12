with import <nixpkgs> {};

pkgs.python35Packages.buildPythonPackage rec {
  name = "PolymerCpp-${version}";
  version = "0.0.1";

  src = ./.;

  propagatedBuildInputs = [ pkgs.python35Packages.tkinter
  			    pkgs.python35Packages.numpy
                            pkgs.python35Packages.matplotlib
  ];

  meta = {
    description = "3D wormlike chain generator for Python and written in C++";
    homepage = "https://github.com/MStefko/PolymerCpp";
  };
}