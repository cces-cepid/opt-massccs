# Optimization parameters for Nitrogen buffer gas

## The repository contents:
- CCS_exp_CO2.dat
- diff_evol_ompc.py
- bounds.dat
- structures in extender xyz format
- log file for each CCS calculation for each molecue. 

### File: CCS_exp_N2.dat

List of molecules and experimental CCS for carbon dioxide buffer gas. The molecules are:
- cocaine.xyz 167.08
- BE.xyz      167.24
- CE.xyz      170.62
- EME.xyz     138.12
- AM.xyz      135.77
- MA.xyz      132.46
- EA.xyz      133.92
- MDA.xyz     145.54
- MDMA.xyz    143.42
- MDEA.xyz    144.06
- GK.xyz      145.43
- AK.xyz      153.82
- VK.xyz      160.60
- GHK.xyz     179.06

### File: bounds.dat

Define the minimal and maximal values of epsilon and sigma parameter for each atom type, for example:

```bash
C 12.011 0.08 0.18 3.0 4.8
N 14.007 0.05 0.18 3.0 4.2
H 1.008  0.01 0.045 1.0 3.0
O 15.999 0.05 0.18 3.0 4.2
```

### File: diff_evol_ompc.py

Python script to perform optimization using differential evolution. Use in this directory the massccs exeutavel and OMPC container. 

```bash
python diff_evol_ompc.py CCS_exp_CO2.dat bounds.dat
```
In this example is used 4 nodes for CCS calculations, see lines 112 and 199, modified this line to perform calculations using different nodes.
