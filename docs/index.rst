.. PolymerCpp documentation master file, created by
   sphinx-quickstart on Wed May 31 07:03:27 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

PolymerCpp
==========

3D wormlike chain generator for Python and written in C++.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

About
=====

PolymerCpp is a small program for generating three-dimensional
wormlike chains (WLC), a common and relatively simple model in polymer
physics. A WLC describes a semi-flexible polymer, i.e. one that is
rigid over short length scales and flexible over long ones. The
characteristic length scale that separates these two regimes is known
as the persistence length.

PolymerCpp provides a number of Python functions that exposes the C++
routines for generating 3D WLCs, including

1. an infinitesimally thin WLC
2. self-avoiding WLCs
3. self-avoiding WLCs with the Rosenbluth method

PolymerCpp was written by `Marcel Stefko
<https://github.com/MStefko>`_ and `Kyle M. Douglass
<https://github.com/kmdouglass>`_ in the `Laboratory of Experimental
Biophysics <http://leb.epfl.ch/>`_ for modeling chromatin
conformations in eukaryotic cells.

User Guide
==========

Installation
------------

Nix
+++

**TODO**
	     
Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
