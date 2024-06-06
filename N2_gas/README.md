# Optimization parameters for Nitrogen buffer gas

## The repository contents:
- CCS_exp_N2.dat
- diff_evol_ompc.py
- bounds.dat
- structures in extender xyz format
- log file for each CCS calculation for each molecue. 

### File: CCS_exp_N2.dat

List of molecules and experimental CCS for nitrogen buffer gas. The molecules are:
- 05_Phenacetin.xyz          140.5
- 07_Alprenolol.xyz          157.5
- 09_Trimethoprim.xyz        171.2
- 10_Ondansetron.xyz         172.7
- 11_Cinchonine_01.xyz       166.3
- 12_Sulpiride_01.xyz        183.3
- 18_Verapamil.xyz           210.0
- 20_Reserpine.xyz           254.3
- acetaminophen.xyz          131.1
- acetylcholine.xyz          127.8
- betamethasone.xyz          189.6
- c60.xyz                    213.1
- c70.xyz                    231.4
- choline.xyz                115.4
- dexamethasone.xyz          190.7
- N-ethylaniline.xyz         124.5
- paracetamol.xyz            131.1
- pyrene.xyz                 135.0
- TMA.xyz                    107.4
- tryphenylene.xyz           143.3
- TtEA.xyz                   122.2

### File: bounds.dat

Define the minimal and maximal values of epsilon and sigma parameter for each atom type, for example:
```bash 
C 12.011 0.06 0.14 3.2 4.0
N 14.007 0.07  0.1  3.8 4.8
H 1.008 0.018 0.05  0.8 2.4
O 15.999 0.05  0.068 3.0 4.0
S 32.06 0.11 0.15 3.0 4.6
P 30.9738 0.14 0.16 3.0 4.6
F 18.998 0.04 0.065 2.8 3.8
```

### File: diff_evol_ompc.py

Python script to perform optimization using differential evolution. Use in this directory the massccs exeutavel and OMPC container. 

```bash
python diff_evol_ompc.py CCS_exp_N2.dat bounds.dat
```
In this example is used 4 nodes for CCS calculations, see lines 112 and 199, modified this line to perform calculations using different nodes.
