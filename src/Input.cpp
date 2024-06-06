/*
 * This program is licensed granted by STATE UNIVERSITY OF CAMPINAS - UNICAMP ("University")
 * for use of MassCCS software ("the Software") through this website
 * https://github.com/cces-cepid/massccs (the "Website").
 *
 * By downloading the Software through the Website, you (the "License") are confirming that you agree
 * that your use of the Software is subject to the academic license terms.
 *
 * For more information about MassCCS please contact: 
 * skaf@unicamp.br (Munir S. Skaf)
 * guido@unicamp.br (Guido Araujo)
 * samuelcm@unicamp.br (Samuel Cajahuaringa)
 */ 

#include "headers/Input.h"

/**
 * Constructor
 * @param argv
 */
Input::Input(char const *argv) {
  readInputFile(argv);
}

/**
 *
 * @param input
 */
void Input::parseFile(char const *input) {
  ifstream inFile;

  inFile.open(input);

  string inpString((istreambuf_iterator<char>(inFile)),
                        istreambuf_iterator<char>());

  const char *json = inpString.c_str();

  d.Parse(json);
  inFile.close();
}

/**
 *
 * @param input
 */
void Input::readInputFile(char const *input) noexcept(false) {
  parseFile(input);

  if (!d.IsObject())
    throw exception();

  /********************************************************
   * Required Parameters
   ********************************************************/
  string option;
  try {
    option = "targetList";
    checkInput(option, 1);
    targetList = d[option.c_str()].GetString();

  } catch (std::invalid_argument &ex) {
    perror(ex.what());
    std::exit(1);
  }

  /********************************************************
   * input parameters
   ********************************************************/
  // number of molecules probes
  if (d.HasMember("numberProbe")) {
    nProbe = d["numberProbe"].GetUint();
  } else {
    nProbe = NPROBE;
  }

  // numbers of CCS calculations
  if (d.HasMember("nIter")) {
    nIter = d["nIter"].GetUint();
  } else {
    nIter = NITER;
  }

  // seed Number to Mersenne Twister - (pseudo)Random number generation
  if (d.HasMember("seed")) {
    seed = d["seed"].GetUint();
  } else {
    seed = SEED;
  }

  // numbers of threads
  if (d.HasMember("nthreads")) {
    nthreads = d["nthreads"].GetUint();
  } else {
    nthreads = omp_get_max_threads();
  }

  // time step in fs
  if (d.HasMember("dt")) {
    dt = d["dt"].GetDouble();
  } else {
    dt = TIMESTEP;
  }
  
  // temperature in Kelvin
  if (d.HasMember("Temp")) {
    temperatureTarget = d["Temp"].GetDouble();
  } else {
    temperatureTarget = TEMPERATURE;
  } 

  // skin of cell size
  if (d.HasMember("skin")) {
    skin = d["skin"].GetDouble();
  } else {
    skin = SKIN;
  }

  // gas buffer type
  if (d.HasMember("GasBuffer")) {
    gas_buffer_str = d["GasBuffer"].GetString();
    if (gas_buffer_str == "He") {
      gas_buffer_flag = 1; 
    } else if (gas_buffer_str == "N2") {
      gas_buffer_flag = 2;
    } else if (gas_buffer_str == "CO2") {
      gas_buffer_flag = 3;  
    } else {
      printf("only available the follow gas buffer types: He N2 CO2\n");
      exit (EXIT_FAILURE);  
    }
  } else {
    gas_buffer_str = "He";
    gas_buffer_flag = 1;
  }

  // equipotential flag
  if (d.HasMember("Equipotential")) {
    equipotential_str = d["Equipotential"].GetString();
    if (equipotential_str == "yes") {
      equipotential_flag =1; 
    } else if (equipotential_str == "no") {
      equipotential_flag = 0;
    } else {
      printf("need to choice Equipotential: yes or no\n");
      exit (EXIT_FAILURE);
    }	    
  } else {
    equipotential_str = "no";
    equipotential_flag = 0;
  }

  // short range flag
  if (d.HasMember("Short-range cutoff")) {
    short_range_str = d["Short-range cutoff"].GetString();
    if (short_range_str == "yes") {
      short_range_cutoff =1; 
    } else if (short_range_str == "no") {
      short_range_cutoff = 0;
    } else {
      printf("need to choice Short-range cutoff: yes or no\n");
      exit (EXIT_FAILURE);
    }	    
  } else {
    short_range_str = "yes";
    short_range_cutoff = 1;
  }

  // lennard-jones cutoff
  if (d.HasMember("LJ-cutoff")) {
    lj_cutoff = d["LJ-cutoff"].GetDouble();
  } else {
    lj_cutoff = SHORT_CUTOFF;
  }
 
  // apply long-range interactions
  if (d.HasMember("Long-range forces")) {
    long_range = d["Long-range forces"].GetString();
    if (long_range == "yes") {
      long_range_flag =1; 
    } else if (long_range == "no") {
      long_range_flag = 0;
    } else {
      printf("need specify only yes or no for coulomb interactions\n");
      exit (EXIT_FAILURE);
    }	    
  } else {
    long_range = "no";
    long_range_flag = 0;
  }

  // long range flag
  if (d.HasMember("Long-range cutoff")) {
    long_range_str = d["Long-range cutoff"].GetString();
    if (long_range_str == "yes") {
      long_range_cutoff =1; 
    } else if (long_range_str == "no") {
      long_range_cutoff = 0;
    } else {
      printf("need to choice Long-range cutoff: yes or no\n");
      exit (EXIT_FAILURE);
    }	    
  } else {
    long_range_str = "yes";
    long_range_cutoff = 1;
  }

  // coulomb cutoff
  if (d.HasMember("Coul-cutoff")) {
    coul_cutoff = d["Coul-cutoff"].GetDouble();
  } else {
    coul_cutoff = LONG_CUTOFF;
  }

  // polarizability flag
  if (d.HasMember("polarizability")) {
    polarizability_str = d["polarizability"].GetString();
    if (polarizability_str == "yes") {
      polarizability_flag =1; 
    } else if (polarizability_str == "no") {
      polarizability_flag = 0;
    } else {
      printf("need to choice polarizability: yes or no\n");
      exit (EXIT_FAILURE);
    }	    
  } else {
    if (gas_buffer_flag == 1) {
      polarizability_str = "yes";
      polarizability_flag = 1;
    } else {
      polarizability_str = "no";
      polarizability_flag = 0;
    }  
  }

  // alpha
  if (d.HasMember("alpha")) {
    alpha = d["alpha"].GetDouble();
  } else {
    if (polarizability_flag == 1) {
      if (gas_buffer_flag == 1) alpha = ALPHA_HE;
      else if (gas_buffer_flag == 2) alpha = ALPHA_N2;
      else if (gas_buffer_flag == 3) alpha = ALPHA_CO2;
      else alpha = 0.0;
    }  
  }

  // force field flag
  if (d.HasMember("force-field")) {
    user_ff = d["force-field"].GetString();
    user_ff_flag = 1;
  } else {
    user_ff_flag = 0;
  }

  // mutation parameter
  if (d.HasMember("mutation")) {
    mut = d["mutation"].GetDouble();
  } else {
    mut = 0.8;
  }

  // crossover parameter 
  if (d.HasMember("crossover")) {
    cross = d["crossover"].GetDouble();
  } else {
    cross = 0.7;
  }

  // numbers of generations for differential evolution
  if (d.HasMember("ngenerations")) {
    ngen = d["ngenerations"].GetUint();
  } else {
    ngen = 100;
  }

  // numbers of individuals
  if (d.HasMember("npopulations")) {
    npop = d["npopulations"].GetUint();
  } else {
    npop = 10;
  }

  // bounds file
  if (d.HasMember("bounds")) {
    bounds= d["bounds"].GetString();
  }

  // numbers of individuals
  if (d.HasMember("ngroups")) {
    ngroups = d["ngroups"].GetUint();
  } else {
    ngroups = 5;
  }

  // simulated annealing temperature parameter 
  if (d.HasMember("Tsa")) {
    Tsa = d["Tsa"].GetDouble();
  } else {
    Tsa = 1000.;
  }

  // simulated annealing cooling rate parameter 
  if (d.HasMember("cooling-rate")) {
    cooling_rate = d["cooling-rate"].GetDouble();
  } else {
    cooling_rate = 0.7;
  }
 

  // read parameters
  /*
  if (d.HasMember("atomicParameters")) {
    const Value& a = d["atomicParameters"];
    for (SizeType i = 0; i < a.Size(); i++) {
      cout << "atom: " << a[i]["atomType"].GetString() << endl;
      cout << "mass: " << a[i]["mass"].GetDouble() << endl;
      cout << "epsilon: " << a[i]["epsilon"][0].GetDouble() << endl;
      cout << "epsilon: " << a[i]["epsilon"][1].GetDouble() << endl;
      cout << "sigma: " << a[i]["sigma"][0].GetDouble() << endl;  
      cout << "sigma: " << a[i]["sigma"][1].GetDouble() << endl;   
    }
  }*/ 

 
}

// 1 = string
// 2 = int
// 3 = double
// 4 = array
/**
 *
 * @param option
 * @param type :
 *
 * * 1 = string
 * * 2 = int
 * * 3 = double
 * * 4 = array
 */
void Input::checkInput(string option, int type) {
  const char *opt = option.c_str();
  if (!d.HasMember(opt))
    throw invalid_argument("Error input.json: " + option + " not found");
  switch (type) {
  case 1:
    if (!d[opt].IsString())
      throw invalid_argument("Error input.json: " + option +
                                  " format error");
    break;
  case 2:
    if (!d[opt].IsInt())
      throw invalid_argument("Error input.json: " + option +
                                  " format error");
    break;
  case 3:
    if (!d[opt].IsDouble())
      throw invalid_argument("Error input.json: " + option +
                                  " format error");
    break;
  case 4:
    if (!d[opt].IsArray())
      throw invalid_argument("Error input.json: " + option +
                                  " format error");
    break;

  default:
    break;
  }
}

/**
 *
 */
void Input::printReadInput() {
  cout << "*********************************************************"
            << endl;
  cout << "INPUT:: Simulation Parameters " << endl;
  cout << "*********************************************************"
            << endl;
  cout << "target list                      : " << targetList << endl;
  cout << "number of probe                  : " << nProbe << endl;
  cout << "number of iterarions             : " << nIter << endl;
  cout << "number of threads                : " << nthreads << endl;
  cout << "seed number                      : " << seed << endl;
  cout << "gas buffer                       : " << gas_buffer_str << endl;
  cout << "Target Temperature (K)           : " << temperatureTarget << endl;
  cout << "timestep (fs)                    : " << dt << endl;
  cout << "Skin cell size (Ang)             : " << skin << endl;
  cout << "Equipotential                    : " << equipotential_str << endl;
  cout << "Cut short-range interaction      : " << short_range_str << endl;
  if (short_range_cutoff == 1) {
  cout << "LJ cutoff (Ang)                  : " << lj_cutoff << endl;
  }
  cout << "Apply long-range interaction     : " << long_range << endl;
  cout << "Cut long-range interaction       : " << long_range_str << endl;
  if (long_range_cutoff == 1) {
  cout << "Coulomb cutoff (Ang)             : " << coul_cutoff << endl;
  }
  cout << "Apply induced-dipole interaction : " << polarizability_str << endl;
  if (polarizability_flag == 1) {
  cout << "alpha (Ang^3)                    : " << alpha << endl;
  }
  if (user_ff_flag == 1) {
    cout << "force-field                      : " << user_ff << endl;
  }
  //cout << "mutation                         : " << mut << endl;
  //cout << "crossover                        : " << cross << endl;
  //cout << "numbers of generations           : " << ngen << endl;
  //cout << "numbers of individuals           : " << npop << endl;
  //cout << "bounds parameters                : " << bounds << endl;
  //cout << "numbers of groups                : " << ngroups << endl;
  //cout << "simulated annealing temperature  : " << Tsa << endl;
  //cout << "cooling rate                     : " << cooling_rate << endl;

}
