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
  
  double end_simulation = omp_get_wtime();

  FILE* ccs_out;
  ccs_out = fopen("ccs.out","w");

  cout << "*********************************************************" << endl;
  cout << "Total time execution " << end_simulation - start_simulation << endl;
  cout << "CCS information " << endl;
  cout << "*********************************************************" << endl;
  cout << "worker molecule mass charge ccs ccs_err time " << endl;

  for (int i = 0; i < nfiles; i++) {
    cout << work_id[i] << "  " << targets[i] << "  " << mass_teor[i] << "  " << charge_teor[i] << "  " << ccs_teor[i] << "  " << ccs_err_teor[i] << "  " << times_teor[i] << endl;
    fprintf(ccs_out,"%s %f\n", targets[i].c_str(), ccs_teor[i]);
  }
  fclose(ccs_out);

  cout << "Program finished..." << endl;
  return 0;

}
