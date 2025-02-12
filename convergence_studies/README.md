## Introduction

This folder contains files for the mesh convergence study of fluid flow simulations.
One dependency for this is [splinepy](https://github.com/tataratat/splinepy).

### Overview of files

- `stokes_studies.json`: case descriptions of analytical solutions using the Stokes
equations
- `export_helpers.py`: splinepy helper functions to export to g+smo-compatible xml-files
- `createXML.py`: creates xml-files for the cases