# PolymerCpp

A C++ implementation of a Worm-like Chain generation algorithm designed to simulate telomere conformation ensembles and calculate their gyration radii. Methods implemented are:
 - Worm-like Chain (WLC) - a Monte Carlo simulation of the Kratky-Porod chain model.
 - Self-avoiding Worm-like Chain (SAWLC) - an adaptation of the WLC model where the generated chains do not overlap with themselves.
 - Rosenbluth SAWLC - adaptation of the SAWLC algorithm which removes bias from the ensemble, but is too slow for practical application - can be used to estimate the influence of this bias.

Parameters of the simulation are set up in the 'main()' function before compiling. The compiler command is:

```bash
g++ -std=c++11 -O2 -I src/Eigen src/*.cpp -o Debug/PolymerCpp -fopenmp
```

No additional dependencies should be required. External libraries used are:
 - [Eigen](http://eigen.tuxfamily.org/) - linear algebra
 - [Stopwatch](https://code.google.com/p/cpp-stopwatch/) - Code profiling and timing

Simulation results are saved in raw text format into specified files, contents of which are completely overwritten.
