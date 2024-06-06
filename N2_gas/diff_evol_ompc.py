#!/usr/bin/env python
import sys,os,string,math
import numpy as np
import json
import time

start = time.time()

# read name of molecule and experimental CCS
if len(sys.argv) != 3:
  print('Not found file list of molecules and its experimental CCS values' %sys.argv[0])
  sys.exit(1)

MoleculeList = str(sys.argv[1])
files = open(MoleculeList, 'r')
lines = files.readlines()

CCS_exp = np.zeros(len(lines)) 
molecules = []
i = 0
for line in lines:
  molecules.append(line.split()[0])
  CCS_exp[i] = float(line.split()[1])
  i += 1

files.close()

#max_bound, min_bound of force field parameters for chemical elements such as: H, C, N, O, F, S, P
elements = []
mass = []
parameter_range = str(sys.argv[2])
pfile = open(parameter_range, 'r')
lines = pfile.readlines()
bounds = [] 

for line in lines:
  elements.append(line.split()[0])
  mass.append(line.split()[1])
  # epsilon bounds
  eps_min = float(line.split()[2]) # min eps
  eps_max = float(line.split()[3]) # max eps
  bounds.append((eps_min, eps_max))
  # sigma bounds
  sig_min = float(line.split()[4]) # min sig
  sig_max = float(line.split()[5]) # max sig  
  bounds.append((sig_min,sig_max))

# differential evolution paramters
mut=0.8
crossp=0.7 
popsize=8
its=400

dimensions = len(bounds)
pop = np.random.rand(popsize, dimensions)
min_b, max_b = np.asarray(bounds).T
diff = np.fabs(min_b - max_b)
pop_denorm = min_b + pop * diff

# massccs parameters
nprobe = 10000
nthreads = 48
niter = 10
seed = 2104
dt = 10.0
Temperature = 298.0
skin = 0.01

CCS_teor = np.zeros(len(CCS_exp)) 
# objetive function values 
fitness = np.zeros(popsize)

os.system('export OMP_NUM_THREADS=48')
os.system('export OMPCLUSTER_SCHEDULER=roundrobin')
os.system('export OMPCLUSTER_BLOCKING_SCHEDULER=1')

# calculate the objetive function for each population
for i in range(popsize):
  # 1. create a force field file 
  ff = open("N2.ff", 'w')
  ff.write("%i\n" % len(elements))
  ff.write("%s\n" % "Symbol  mass (amu)  epislon (kcal/mol)  sigma(Angstroms)")

  for j in range(len(elements)):
    ff.writelines("%s %s %f %f\n" % (elements[j], mass[j], pop_denorm[i][2*j], pop_denorm[i][2*j+1]))
  ff.close()
  
  # CCS calculation using OMPC
  input = {
    "targetList": "CCS_exp_N2.dat",
    "force-field" : "N2.ff",
    "numberProbe" : nprobe,
    "nthreads" : nthreads,
    "nIter" : niter,
    "seed" : seed,
    "dt" : dt,
    "Temp" : Temperature,
    "skin" : skin,
    "GasBuffer" : "N2",
    "Equipotential" : "yes",
    "Short-range cutoff" : "no",
    "LJ-cutoff" : 12.0,
    "Long-range forces" : "yes",
    "Long-range cutoff" : "no",
    "Coul-cutoff" : 25.0,
    "polarizability" : "yes"
  }

  with open("input.json", "w") as outfile:
    json.dump(input, outfile)
  # run massccs-opmc using 4 nodes 
  os.system('mpirun --np 4 apptainer exec --nv runtime_latest.sif ./massccs input.json > log.massccs 2>&1')

  # copy the CCS_teor
  ccs = open("ccs.out", 'r')
  lines = ccs.readlines()
  j = 0
  for line in lines:
    CCS_teor[j] = float(line.split()[1])
    j += 1;
  ccs.close() 
  
  # calculation objetive function
  fitness[i] = sum((1-CCS_teor/CCS_exp)**2)
  
best_idx = np.argmin(fitness)
best = pop_denorm[best_idx]

# print evolution of objetive function
evol = open("fobj.dat", 'w')
evol.write("%s\n" % "#differentail evolution iteration")
evol.writelines("%i %g " % (0, fitness[best_idx]))
for j in range(len(elements)):
  evol.writelines("%f %f " % (best[2*j], best[2*j+1]))
evol.writelines("\n")
evol.flush()

fobj = []
pop_evol = []
fobj.append(fitness[best_idx])
pop_evol.append(best)

# print the population 
generation = open("population.dat", 'w')
generation.write("%s\n" % "#iteration population")
generation.writelines("%i " % (0))
for i in range(popsize):
   for j in range(len(elements)): 
     generation.writelines("%f %f " % (pop_denorm[i][2*j], pop_denorm[i][2*j+1]))

generation.writelines("\n")
generation.flush()

# differential evolution run
for i in range(its):
  for j in range(popsize):
    
    idxs = [idx for idx in range(popsize) if idx != j]
    a, b, c = pop[np.random.choice(idxs, 3, replace = False)]
    mutant = np.clip(a + mut * (b - c), 0, 1)
    cross_points = np.random.rand(dimensions) < crossp
    if not np.any(cross_points):
      cross_points[np.random.randint(0, dimensions)] = True
    trial = np.where(cross_points, mutant, pop[j])
    trial_denorm = min_b + trial * diff
    
    # calculate the CCS for trial_denorm 
    ff = open("N2.ff", 'w')
    ff.write("%i\n" % len(elements))
    ff.write("%s\n" % "Symbol  mass (amu)  epislon (kcal/mol)  sigma(Angstroms)")
    for k in range(len(elements)):    
      ff.writelines("%s %s %f %f\n" % (elements[k], mass[k], trial_denorm[2*k], trial_denorm[2*k+1]))
    ff.close()

    # CCS calculation using OMPC
    input = {
      "targetList": "CCS_exp_N2.dat",
      "force-field" : "N2.ff",
      "numberProbe" : nprobe,
      "nthreads" : nthreads,
      "nIter" : niter,
      "seed" : seed,
      "dt" : dt,
      "Temp" : Temperature,
      "skin" : skin,
      "GasBuffer" : "N2",
      "Equipotential" : "yes",
      "Short-range cutoff" : "no",
      "LJ-cutoff" : 12.0,
      "Long-range forces" : "yes",
      "Long-range cutoff" : "no",
      "Coul-cutoff" : 25.0,
      "polarizability" : "yes"  
    }

    with open("input.json", "w") as outfile:
      json.dump(input, outfile)
    # run massccs-opmc using 4 nodes 
    os.system('mpirun --np 4 apptainer exec --nv runtime_latest.sif ./massccs input.json > log.massccs 2>&1')

    # copy the CCS_teor
    ccs = open("ccs.out", 'r')
    lines = ccs.readlines()
    k = 0
    for line in lines:
      CCS_teor[k] = float(line.split()[1])
      k += 1;
    ccs.close()
   
    # calculation objetive function
    f = sum((1-CCS_teor/CCS_exp)**2)

    if f < fitness[j]:
      fitness[j] = f
      pop[j] = trial
      if f < fitness[best_idx]:
        best_idx = j
        best = trial_denorm 

  fobj.append(fitness[best_idx])
  pop_evol.append(best)

  evol.writelines("%i %g " % (i+1, fitness[best_idx]))
  for jj in range(len(elements)):
    evol.writelines("%f %f " % (best[2*jj], best[2*jj+1]))
  evol.writelines("\n")
  evol.flush()

  pop_denorm = min_b + pop * diff

  generation.writelines("%i " % (i+1))
  for ii in range(popsize):
    for jj in range(len(elements)):
      generation.writelines("%f %f " % (pop_denorm[ii][2*jj], pop_denorm[ii][2*jj+1]))
  generation.writelines("\n")
  generation.flush()

evol.close()
generation.close()

end = time.time()
print('time (s) = ', end - start)
