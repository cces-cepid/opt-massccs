/*
 * This program is licensed granted by STATE UNIVERSITY OF CAMPINAS - UNICAMP ("University")
 * for use of MassCCS software ("the Software") through this website
 * https://github.com/cces-cepid/MassCCS (the "Website").
 *
 * By downloading the Software through the Website, you (the "License") are confirming that you agree
 * that your use of the Software is subject to the academic license terms.
 *
 * For more information about MassCCS please contact: 
 * skaf@unicamp.br (Munir S. Skaf)
 * guido@unicamp.br (Guido Araujo)
 * samuelcm@unicamp.br (Samuel Cajahuaringa)
 */

#include "headers/System.h"
#include "headers/Input.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <fstream>
#include <stdio.h>
#include <cstring>
#include <unistd.h>
#include "omp.h"
#include <cstdlib>

using namespace std;

int nfiles;
string *targets;
double *ccs_exp;
string *elements;
double *mass;
double *eps_max;
double *eps_min;
double *sig_max;
double *sig_min;
int dimension;

void readlist(string &targetList) {
  string sstring, recordName;
  ifstream infile;
  string type, line;
  stringstream ss;

  /* Check the number of target files */
  nfiles = 0;
  infile.open (targetList);
  if (infile.is_open()){
    while(getline(infile,sstring)) // To get you all the lines.
    {
      istringstream sstream(sstring);
      sstream >> recordName;
      ++nfiles;
    }
    infile.close();
  }

  //cout << "nfiles: " << nfiles << endl;
  double ccs;
  targets = new string[nfiles];
  ccs_exp = new double[nfiles];

  int i = 0;
  infile.open(targetList);
  if (infile.is_open()) {
    while(getline(infile,sstring))
    {
      istringstream sstream(sstring);
      sstream >> recordName >> ccs;
      targets[i] = recordName;
      ccs_exp[i] = ccs;
      //cout << "targets: " << targets[i] << " ccs: " << ccs_exp[i] << endl;
      ++i;
    }
  }
}

void readbounds(string &bounds) {
  string sstring, recordName;
  ifstream infile;
  string type, line;
  stringstream ss;

  /* Check the number of target files */
  dimension = 0;
  infile.open (bounds);
  if (infile.is_open()){
    while(getline(infile,sstring)) // To get you all the lines.
    {
      istringstream sstream(sstring);
      sstream >> recordName;
      ++dimension;
    }
    infile.close();
  }

  elements = new string[dimension];
  mass = new double[dimension]; 
  eps_max = new double[dimension];
  eps_min = new double[dimension];
  sig_max = new double[dimension];
  sig_min = new double[dimension];
   
  
  double e_max, e_min, s_max, s_min, mi;
  int i = 0;
  infile.open(bounds);
  if (infile.is_open()) {
    while(getline(infile,sstring))
    {
      istringstream sstream(sstring);
      sstream >> recordName >> mi >> e_min >> e_max >> s_min >> s_max;
      elements[i] = recordName;
      mass[i] = mi;
      eps_min[i] = e_min;
      eps_max[i] = e_max;
      sig_min[i] = s_min;
      sig_max[i] = s_max;
      ++i;
    }
  }
  //for (int i = 0; i < dimension; i++) cout << elements[i] << " " << mass[i] << " " << eps_min[i] << " " << eps_max[i] << " " << sig_min[i] << " " << sig_max[i] << endl;
}

void write_ff(double individual[], string &user_ff, int dimension) {
  ofstream outfile;
   
  string def_elements[7] = {"C", "N", "H", "O", "S", "P", "F"};
  double def_mass[7] = {12.011, 14.007, 1.008, 15.999, 32.06, 30.9738, 18.9984};
  double def_eps[7] = {0.0824736, 0.0758527, 0.0362711, 0.062323, 0.138032, 0.145372, 0.04649};
  double def_sig[7] = {3.2255, 3.5719, 1.8986, 3.075, 3.4237, 3.47, 3.1285};

  int k = 0;
  //int nn = elements.size();
  
  for (int i = 0; i < dimension/2; i++) {
    string chem = elements[i];
    if (chem == "C") {
      def_eps[0] = individual[0];//2*k];
      def_sig[0] = individual[1];//2*k+1];
      //++k;
    } else if (chem == "N") {
      def_eps[1] = individual[2];//2*k];
      def_sig[1] = individual[3];//2*k+1];
      //++k;
    } else if (chem == "H") {
      def_eps[2] = individual[4];//2*k];
      def_sig[2] = individual[5];//2*k+1];
      //++k;
    } else if (chem == "O") {
      def_eps[3] = individual[6];//2*k];
      def_sig[3] = individual[7];//2*k+1];
      //++k;       
    } else if (chem == "S") {
      def_eps[4] = individual[8];//2*k];
      def_sig[4] = individual[9];//2*k+1];
      //++k;
    } else if (chem  == "P") {
      def_eps[5] = individual[10];//2*k];
      def_sig[5] = individual[11];//2*k+1];
      //++k;
    } else if (chem == "F") {
      def_eps[6] = individual[12];//2*k];
      def_sig[6] = individual[13];//2*k+1];
      //++k;
    }
  } 

  outfile.open(user_ff);
  outfile << 7 << endl;
  outfile << "Symbol  mass (amu)  epislon (kcal/mol)  sigma(Angstroms)" << endl;
  for (int i = 0; i < 7; i++) outfile << def_elements[i] << " " << def_mass[i] << " " << def_eps[i] << " " << def_sig[i] << endl;
  outfile.close(); 
}
 
double fitness(double ccs_exp[], double ccs_teor[], int nfiles) {
  double f = 0.0;
  for (int i = 0; i < nfiles; i++) f += pow(1.0 - ccs_teor[i]/ccs_exp[i],2);
  return f;
}

int min_idx(double fobj[], int npop) {
  int idx;
  double fmin = 10000000.0;

  for (int i = 0; i < npop; i++) {
    if (fobj[i] < fmin) {
      fmin = fobj[i];
      idx = i;
    }
  } 
  return idx; 
}

void select(int &a, int &b, int &c, int it, int npop) {
  int idx_1[npop-1];
  int k = 0;

  for (int i = 0; i < npop; i++) {
    if (i != it) {
      idx_1[k] = i;
      ++k;
    } 
  }
   
  // pick random number
  int i1 = rand() % (npop-1);
  a = idx_1[i1];

  int idx_2[npop-2];
  k = 0;

  for (int i = 0; i < npop-1; i++) {
    if (idx_1[i] != a) { 
      idx_2[k] = idx_1[i];
      ++k;
    }
  }
   
  // pick random number
  int i2 = rand() % (npop-2);
  b = idx_2[i2];

  int idx_3[npop-3];
  k = 0;

  for (int i = 0; i < npop-2; i++) {
    if (idx_2[i] != b) {
      idx_3[k] = idx_2[i];
      ++k;
    }
  }
  
  // pick random number
  int i3 = rand() % (npop-3);
  c = idx_3[i3];
  
}


int main(int argc, char *argv[]) {
  // check if the input is right
  if (argc < 2) {
    cout << "Usage: " << argv[0] << " input.json" << endl;
    exit(1);
  }
   
  unsigned int seed, nProbe, nIter, equipotential_flag, gas_buffer_flag, nthreads;
  unsigned int short_range_cutoff, long_range_flag, long_range_cutoff, polarizability_flag, user_ff_flag;
  double temperatureTarget, dt, skin, lj_cutoff, coul_cutoff, alpha, mut, cross;
  unsigned int force_type, ngen, npop;
  string targetList, user_ff, bounds;
  /**
   * Read the input file
   */
  double start_simulation = omp_get_wtime(); 
 
  Input *input;
  input = new Input(argv[1]);

  input->printReadInput();

  // initialize variables from input values
  nProbe = input->nProbe;                           // numbers of gas buffers 
  nIter = input->nIter;                             // numbers of CCS calculations
  seed = input->seed;                               // random number seed 
  nthreads = input->nthreads;                       // number of threads
  targetList = input->targetList;                   // xyz or pqr file of molecule target 
  dt = input->dt;                                   // time step in fs 
  temperatureTarget = input->temperatureTarget;     // temperature in Kelvin
  gas_buffer_flag = input->gas_buffer_flag;         // He = 1, N2 = 2
  equipotential_flag = input->equipotential_flag;
  skin = input->skin;                               // skin cell-size   
  short_range_cutoff = input->short_range_cutoff;   // yes = 1 and not = 0 for cut lennard-jones interacion
  lj_cutoff = input->lj_cutoff;                     // lennard-jones cutoff   
  long_range_flag =  input->long_range_flag;        // apply coulomb interaction 
  long_range_cutoff = input->long_range_cutoff;     // yes = 1 and not = 0 for cut coulomb interacion
  coul_cutoff = input->coul_cutoff;                 // coulomb cutoff
  polarizability_flag = input->polarizability_flag; // apply induced-dipole interaction
  alpha = input->alpha;                             // polarizability
  user_ff_flag = input->user_ff_flag;               // user force field
  user_ff = input->user_ff;                         // name of force field 
  mut = input->mut;
  cross = input->cross;
  ngen = input->ngen;
  npop = input->npop;
  bounds = input->bounds;

  if (gas_buffer_flag == 1) {
    if (short_range_cutoff == 0 && long_range_flag == 0) {
      // only lennard-jones interaction 
      force_type = 1;
    } else if (short_range_cutoff == 1 && long_range_flag == 0) {
      // only lennard-jones interaction with cutoff
      force_type = 2;
    } else if (short_range_cutoff == 0 && polarizability_flag == 1 && long_range_cutoff == 0) {
      // lennard-jones and induced-dipole interactions
      force_type = 3;
    } else if (short_range_cutoff == 1 && polarizability_flag == 1 && long_range_cutoff == 1) {
      // lennard-jones and induced-dipole interactions with cutoff
      force_type = 4;
    } else {
      force_type = 2;
    }
  } else if (gas_buffer_flag == 2 || gas_buffer_flag == 3){
    if (short_range_cutoff == 0 && long_range_flag == 0 && polarizability_flag == 0) {
    // only lennard-jones interaction 
      force_type = 1;
    } else if (short_range_cutoff == 1 && long_range_flag == 0 && polarizability_flag == 0) {
      // only lennard-jones interaction with cutoff
      force_type = 2;
    } else if (short_range_cutoff == 0 && long_range_flag == 1 && long_range_cutoff == 0 && polarizability_flag == 0) {
      // lennard-jones and coulomb interactions
      force_type = 3;
    } else if (short_range_cutoff == 1 && long_range_flag == 1 && long_range_cutoff == 1 && polarizability_flag == 0) {
      // lennard-jones and coulomb interactions with cutoff
      force_type = 4;
    } else if (short_range_cutoff == 0 && long_range_flag == 1 && long_range_cutoff == 0 && polarizability_flag == 1) {
      // lennard-jones, coulomb and induced-dipole interactions 
      force_type = 5;
    } else if (short_range_cutoff == 1 && long_range_flag == 1 && long_range_cutoff == 1 && polarizability_flag == 1) {
      // lennard-jones, coulomb and induced-dipole interactions with cutoff
      force_type = 6;
    } else {
      force_type = 2;
    }
  } 
 
  // read the list of molecule target
  readlist(targetList);
  
  cout << "targets files " << endl;
  cout << "nfiles: " << nfiles << endl;
  for (int i = 0; i < nfiles; i++) cout << targets[i] << endl;
  
  double ccs_teor[nfiles];
  double ccs_err_teor[nfiles];
  double mass_teor[nfiles];
  double charge_teor[nfiles];
  double times_teor[nfiles];
  int work_id[nfiles];
   
  for (int i = 0; i < nfiles; i++) {
    ccs_teor[i] = 0.;
    ccs_err_teor[i] = 0.;
    mass_teor[i] = 0.;
    charge_teor[i] = 0.;
    times_teor[i] = 0.;
    work_id[i] = 0;
  }

  readbounds(bounds);
  // populations
  
  dimension *= 2;
  double pop[npop][dimension];
  double pop_denorm[npop][dimension];
  double individual[dimension];

  RandomNumber *rn{};
  rn = new RandomNumber(seed+dimension);
   
  for (int i = 0; i < npop; i++) {
    for (int j = 0; j < dimension; j++) {
        pop[i][j] = rn->getRandomNumber();
    }
  }

  double min_b[dimension], max_b[dimension];   
  for (int j = 0; j < dimension/2; j++) {
    min_b[2*j] = eps_min[j];
    min_b[2*j+1] = sig_min[j];
    max_b[2*j] = eps_max[j];
    max_b[2*j+1] = sig_max[j];
    //cout << min_b[2*j] << " " << max_b[2*j] << endl;
    //cout << min_b[2*j+1] << " " << max_b[2*j+1] << endl;
  }
  
  for (int i = 0; i < npop; i++) {
    for (int j = 0; j < dimension; j++) {
      pop_denorm[i][j] = min_b[j] + abs(max_b[j] - min_b[j])*pop[i][j];
      //cout << "pop denorm: " << pop_denorm[i][j] << endl;
    }
  }

  // convert list of file to vector char
  int nsize = 0;
  int nlocal[nfiles];

  for (int i = 0; i < nfiles; i++) {
    nlocal[i] = targets[i].size();
    nsize += nlocal[i];
  }

  char targets_word[nsize];
  int jj = 0;
  for (int i = 0; i < nfiles; i++) {
    for (int j = 0; j < targets[i].size(); j++) {
      targets_word[jj] = targets[i][j];
      jj += 1;
    }
  }

  int nbegin[nfiles];
  nbegin[0] = 0;

  for (int i = 1; i < nfiles; i++) nbegin[i] = nbegin[i-1] + nlocal[i-1];

  // convert name force field to vector char
  int nff = user_ff.size();
  char user_ff_word[nff];

  for (int i = 0; i < user_ff.size(); i++) user_ff_word[i] = user_ff[i];
  
  int fake_var;

  double fobj[npop];
  
  for (int it = 0; it < npop; it++)  {
    // user force field for each individual
    for (int id = 0; id < dimension; id++) individual[id] = pop_denorm[it][id];

    write_ff(individual, user_ff, dimension);  
       
    #pragma  omp parallel
    #pragma omp single
    {
      for (int i = 0; i < nfiles; i++) {
      #pragma omp target nowait \
       map(to: nProbe, nIter, seed, nthreads, dt, temperatureTarget, gas_buffer_flag, equipotential_flag, skin) \
       map(to: short_range_cutoff, lj_cutoff, long_range_flag, long_range_cutoff, coul_cutoff, polarizability_flag, alpha, user_ff_flag, force_type) \
       map(to: user_ff_word[:nff]) \
       map(to: targets_word[nbegin[i]:nlocal[i]]) \
       map(tofrom: ccs_teor[i:1], ccs_err_teor[i:1], mass_teor[i:1], charge_teor[i:1], times_teor[i:1], work_id[i:1]) \
       depend(in: fake_var)
       {
         string target;
         for (int j = 0; j < nlocal[i]; j++) {
           int k = nbegin[i] + j;
           target.push_back(targets_word[k]);
         }

         string file_ff;
         for (int j = 0; j < nff; j++) {
           file_ff.push_back(user_ff_word[j]);
         }
         work_id[i] = getpid();
   
         size_t lastindex = target.find_last_of(".");
         string name = target.substr(0,lastindex);
         string logname;
         logname = name+".log";   

         double start_conformation = omp_get_wtime();     
         System *system; 
         system = new System(target, logname, nProbe, nIter, seed, nthreads, dt, temperatureTarget, gas_buffer_flag, equipotential_flag, skin, short_range_cutoff, lj_cutoff, 
         long_range_flag, long_range_cutoff, coul_cutoff, polarizability_flag, alpha, file_ff, user_ff_flag, force_type, ccs_teor[i], ccs_err_teor[i], mass_teor[i], charge_teor[i]);
         delete system;
         double end_conformation = omp_get_wtime();
        
         times_teor[i] = end_conformation - start_conformation;
       }
      } 
    }
    // for (int i = 0; i < nfiles; i++) cout << ccs_teor[i] << "  " << ccs_exp[i] << endl;
    fobj[it] = fitness(ccs_exp, ccs_teor, nfiles);
  }
  
  // best idx and best individual
  int best_idx;
  double best[dimension];

  best_idx = min_idx(fobj, npop);

  for (int id = 0; id < dimension; id++) best[id] = pop_denorm[best_idx][id];

  // write fob.dat and generation.dat
  ofstream fobjetive;
  fobjetive.open ("fobj.dat");
  fobjetive << "#differentail evolution iteration" << endl;
  fobjetive << "0 " << fobj[best_idx] << " ";

  for (int id = 0; id < dimension; id++) fobjetive << best[id] << " ";
  fobjetive << "\n";
  fobjetive.flush();
 
  ofstream fgeneration;
  fgeneration.open ("generation.dat");
  fgeneration << "#iteration population" << endl;
  fgeneration << "0 ";

  for (int ii = 0; ii < npop; ii++) {
    for (int id = 0; id < dimension; id++) {
      fgeneration << pop_denorm[ii][id] << " ";
    }
  }
  fgeneration << "\n";
  fgeneration.flush();

   
  int i1, i2, i3;
  double mutation[dimension], trial[dimension], trial_denorm[dimension];
  bool crossover[dimension];
  double fmin;
  double cross_point;
  bool any; 

  // differential evolution
  for (int ii = 0; ii < ngen; ii++)  {
    for (int it = 0; it < npop; it++)  {
      // select 3 random numbers exept the it
      select(i1, i2, i3, it, npop);
    
      // mutation
      for (int id = 0; id < dimension; id++) {
        mutation[id] = pop[i1][id] + mut * (pop[i2][id] - pop[i3][id]);
        if (mutation[id] < 0.) {
          mutation[id] = 0.;
        } else if (mutation[id] > 1.) {
          mutation[id] = 1.;
        }
      }   
   
      // crossover
      for (int id = 0; id < dimension; id++) {
        cross_point = rn->getRandomNumber();
        if (cross_point < cross) {
          crossover[id] = true;
        } else {
          crossover[id] = false;
        }
      }
    
      any = false;
      for (int id = 0; id < dimension; id++) {
        if (crossover[id]) any = true;
      } 
    
      if (!any) crossover[rand() % dimension] = true;
 
      // trial and trial_denorm
      for (int id = 0; id < dimension; id++) {
        if (crossover[id]) trial[id] = mutation[id];
        else trial[id] = pop[it][id];
      }
  
      for (int id = 0; id < dimension; id++) trial_denorm[id] = min_b[id] + abs(max_b[id] - min_b[id])*trial[id];
    
      // write ff
      write_ff(trial_denorm, user_ff, dimension);
  
      // ccs calculations
      #pragma  omp parallel
      #pragma omp single
      {
        for (int i = 0; i < nfiles; i++) {
        #pragma omp target nowait \
         map(to: nProbe, nIter, seed, nthreads, dt, temperatureTarget, gas_buffer_flag, equipotential_flag, skin) \
         map(to: short_range_cutoff, lj_cutoff, long_range_flag, long_range_cutoff, coul_cutoff, polarizability_flag, alpha, user_ff_flag, force_type) \
         map(to: user_ff_word[:nff]) \
         map(to: targets_word[nbegin[i]:nlocal[i]]) \
         map(tofrom: ccs_teor[i:1], ccs_err_teor[i:1], mass_teor[i:1], charge_teor[i:1], times_teor[i:1], work_id[i:1]) \
         depend(in: fake_var)
         {
           string target;
           for (int j = 0; j < nlocal[i]; j++) {
             int k = nbegin[i] + j;
             target.push_back(targets_word[k]);
           }

           string file_ff;
           for (int j = 0; j < nff; j++) {
             file_ff.push_back(user_ff_word[j]);
           }

           work_id[i] = getpid();

           size_t lastindex = target.find_last_of(".");
           string name = target.substr(0,lastindex);
           string logname;
           logname = name+".log";

           double start_conformation = omp_get_wtime();
           System *system;
           system = new System(target, logname, nProbe, nIter, seed, nthreads, dt, temperatureTarget, gas_buffer_flag, equipotential_flag, skin, short_range_cutoff, lj_cutoff,
           long_range_flag, long_range_cutoff, coul_cutoff, polarizability_flag, alpha, file_ff, user_ff_flag, force_type, ccs_teor[i], ccs_err_teor[i], mass_teor[i], charge_teor[i]);
           delete system;
           double end_conformation = omp_get_wtime();

           times_teor[i] = end_conformation - start_conformation;
         }
        }
      }

      fmin = fitness(ccs_exp, ccs_teor, nfiles);

      if (fmin < fobj[it]) {
        fobj[it] = fmin;
        for (int id = 0; id < dimension; id++) pop[it][id] = trial[id];
      
        if (fmin < fobj[best_idx]) {
          best_idx = it;       
          for (int id = 0; id < dimension; id++) best[id] = trial_denorm[id];
        }
      } 
    }

    // k random number between 2 and \sqrt(npop)
    k = 2 + rand() % (int(sqrt(npop)) - 2 + 1)
   
    // caluculate the k clusters center positions and the new populations
    



    //   for (int id = 0; id < dimension; id++) best[id] = pop_denorm[best_idx][id];
    // write fob.dat and generation.dat
    fobjetive << ii+1 << " " << fobj[best_idx] << " ";
    for (int id = 0; id < dimension; id++) fobjetive << min_b[id] + abs(max_b[id] - min_b[id])*pop[best_idx][id] << " ";
    fobjetive << "\n";
    fobjetive.flush();     

    fgeneration << ii+1 << " ";
    for (int ip = 0; ip < npop; ip++) {
      for (int jp = 0; jp < dimension; jp++) {
        pop_denorm[ip][jp] = min_b[jp] + abs(max_b[jp] - min_b[jp])*pop[ip][jp];      
        fgeneration << pop_denorm[ip][jp] << " ";
      }
    }
    fgeneration << "\n";
    fgeneration.flush();
  }

  fobjetive.close();
  fgeneration.close();
  double end_simulation = omp_get_wtime();
  cout << "Total time: " << end_simulation - start_simulation << endl;
  cout << "Program finished..." << endl;
  return 0;
}
