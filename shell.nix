with import <nixpkgs> {};

( let

    PolymerCpp = pkgs.python35Packages.buildPythonPackage rec {
      name = "PolymerCpp-${version}";
      version = "0.0.0";

      src = ./.;

      propagatedBuildInputs = [ pkgs.python35Packages.numpy ];
      
    };

  in pkgs.python35.withPackages (ps: [ ps.numpy
     				       ps.sphinx
				       pkgs.ncurses
				       pkgs.netbeans
				       PolymerCpp])
).env