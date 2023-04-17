#include "simulation_support.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

// Neutrons per gram of Californium 252 released per second
#define CF252_NEUTRONS_PER_GRAM_PER_SEC 23e12
// Number of nanoseconds in a second
#define NS_AS_SEC 1e-9
// Speed of light in m/s
#define C 3e8
// Constant of proportionality
#define R0 1.2e-15
// PI
#define PI 3.14159265359
// Planck's constant multiplied by speed of light
#define HC 1239.84
// Size of a barn in meters squared
#define BARN 1e-28
// The Avogadro constant
#define AVOGADRO 6.022e23
// Amount of energy in MeV released per fission
#define MeV_PER_FISSION 200
// Number of Joules of energy corresponding to 1 MeV
#define JOULES_PER_MEV 1.6e-13
// The scattering cross section of hydrogen
#define SCATTERING_CROSS_SECTION_HYDROGEN 20
// The scattering cross section of oxygen
#define SCATTERING_CROSS_SECTION_OXYGEN 3.8
// The scattering cross section of carbon
#define SCATTERING_CROSS_SECTION_CARBON 4.9
// Density of water in grams per cm^3
#define DENSITY_WATER 1
// Density of deuterium in grams per cm^3
#define DENSITY_DEUTERIUM 1.8
// Density of graphite in grams per cm^3
#define DENSITY_GRAPHITE 2.25
// Probability that water will absorb neutron
#define H2O_ABSORBTION_PROB 0.25
// Probability that deuterium will absorb neutron
#define D2O_ABSORPTION_PROB 0.05
// Probability that graphite will absorb neutron
#define C_ABSORPTION_PROB 0.01
// Mass in kg of one unit of atomic rest mass
#define MASS_ONE_UNIT (1.6726*1e-27)
// Newton meters per electron volt
#define NEWTON_METER_PER_eV (1.60217733*1e-19)

static bool determineNeutronModeratorScattering(int, enum moderator_type_enum, int);
static double getScatteringCrossSection(enum chemical_type_enum);
static double getScatteringCrossSectionOfWater();
static double getScatteringCrossSectionOfDeuterium();
static double getScatteringCrossSectionOfGraphite();
static void generateNeutronVelocityAndEnergy(struct neutron_struct*);
static bool determineNeutronAbsorbtionByFuel(double, struct channel_struct*, enum chemical_type_enum, int, int);
static double generateDecimalRandomNumber();
static double getRandomNumberGeneratorResolution();
static double getNumberAtomsInFuelPellet(struct channel_struct*, enum chemical_type_enum, int);
static double getCrossSectionInBarns(double, enum chemical_type_enum);
static int getAtomicRestMass(enum chemical_type_enum);
static double calculateDeBroglieWavelength(double, int);

/**
 * Given the weight of a sample of Californium 252, will calculate how many
 * neutrons this has produced in a number of nanoseconds
 **/
unsigned long int getNumberNeutronsFromGenerator(double weight, int dt) {
  return CF252_NEUTRONS_PER_GRAM_PER_SEC * weight * (dt * NS_AS_SEC);
}

/**
 * Converts energy in MeV with a specific atomic rest mass to velocity
 * in meters per second. A neutron always has a rest mass of one
 **/
double MeVToVelocity(double MeV, int rest_mass) {
  double mass=rest_mass*MASS_ONE_UNIT;
  double ke=MeV*1e6*NEWTON_METER_PER_eV;
  double vvc=sqrt(1-1/pow(ke/(mass*C*C)-(-1), 2));
  return vvc * C;
}

/**
 * Will determine and handle if a neutron has collided with the moderator, returning true
 * if the neutron is no longer involved in the simulation (i.e. it has been absorbed) or
 * false if it is still present. On scattering the energy of the neutron is updated.
 * Note that for the absorption we use predefined probabilities, that are fairly close to those
 * measured experimentally rather than computing the absorption cross section directly.
 **/
bool determineAndHandleIfNeutronModeratorCollision(struct neutron_struct * neutron, int moderator_weight,
    enum moderator_type_enum moderator_type, int size_z) {
  if (determineNeutronModeratorScattering(moderator_weight, moderator_type, size_z)) {
    double absorption_prob;
    int atomic_mass;
    if (moderator_type == WATER) {
      absorption_prob=H2O_ABSORBTION_PROB;
      atomic_mass=1;
    }
    if (moderator_type == DEUTERIUM) {
      absorption_prob=D2O_ABSORPTION_PROB;
      atomic_mass=2;
    }

    if (moderator_type == GRAPHITE) {
      absorption_prob=C_ABSORPTION_PROB;
      atomic_mass=12;
    }
    if (generateDecimalRandomNumber() >= absorption_prob) {
      // Use 1 here for the mass as its the hydrogen in the H20 that does the slowing
      double mean_log_energy_reduction=2/(1+(atomic_mass/3));
      neutron->energy=neutron->energy/exp(mean_log_energy_reduction);
      return false;
    } else {
      // Neutron absorbed by water
      return true;
    }
  }
  return false;
}

/**
 * Determines and handles if a neutron has collided with an atom of fuel, returning true if
 * so (the neutron will be deactivated from the simulation) and false otherwise. This handles both
 * the U235 and Pu239 fuels (it is assumed that no other fuels fission).
 **/
bool determineAndHandleIfNeutronFuelCollision(double MeV, struct channel_struct * reactor_channel,
    int fuel_pellet, int collision_prob_multiplyer) {
  double deBroglieWavelength=calculateDeBroglieWavelength(MeV, 1);
  // Calculate the cross section in Barns for each chemical and determine if a collision occured
  // if so add the neutron to the atom it collided with.
  if (reactor_channel->contents.fuel_assembly.quantities[fuel_pellet][U235] > 0) {
    if (determineNeutronAbsorbtionByFuel(deBroglieWavelength, reactor_channel, U235, fuel_pellet, collision_prob_multiplyer)) {
      reactor_channel->contents.fuel_assembly.quantities[fuel_pellet][U235]--;
      reactor_channel->contents.fuel_assembly.quantities[fuel_pellet][U236]++;
      return true;
    }
  }
  if (reactor_channel->contents.fuel_assembly.quantities[fuel_pellet][Pu239] > 0) {
   if (determineNeutronAbsorbtionByFuel(deBroglieWavelength, reactor_channel, Pu239, fuel_pellet, collision_prob_multiplyer)) {
      reactor_channel->contents.fuel_assembly.quantities[fuel_pellet][Pu239]--;
      reactor_channel->contents.fuel_assembly.quantities[fuel_pellet][Pu240]++;
      return true;
    }
  }
  return false;
}

/**
 * Handles the fission of Pu240, which will either split into Xenon and
 * Zirconium, releasing 3 neutrons, or go back to Pu239 releasing
 * one neutron.
 **/
int fissionPu240(struct channel_struct * channel, int fuel_pellet) {
  double r=generateDecimalRandomNumber();
  if (r <0.73) {
    channel->contents.fuel_assembly.quantities[fuel_pellet][Xe134]++;
    channel->contents.fuel_assembly.quantities[fuel_pellet][Zr103]++;
    channel->contents.fuel_assembly.quantities[fuel_pellet][Pu240]--;
    return 3;
  } else {
    channel->contents.fuel_assembly.quantities[fuel_pellet][Pu239]++;
    channel->contents.fuel_assembly.quantities[fuel_pellet][Pu240]--;
    return 1;
  }
}

/**
 * Handles the fission of U236, which will either split into BariumgetNumberAtomsInFuelPellet
 * and Krypton releasing 3 neutrons, or Xenon and Strontium releasing
 * 2 neutrons.
 **/
int fissionU236(struct channel_struct * channel, int fuel_pellet) {
  double r=generateDecimalRandomNumber();
  if (r <0.85) {
    channel->contents.fuel_assembly.quantities[fuel_pellet][Ba141]++;
    channel->contents.fuel_assembly.quantities[fuel_pellet][Kr92]++;
    channel->contents.fuel_assembly.quantities[fuel_pellet][U236]--;
    return 3;
  } else {
    channel->contents.fuel_assembly.quantities[fuel_pellet][Xe140]++;
    channel->contents.fuel_assembly.quantities[fuel_pellet][Sr94]++;
    channel->contents.fuel_assembly.quantities[fuel_pellet][U236]--;
    return 2;
  }
}

/**
 * Initialises a neutron based on it's location in the reactor (the channel
 * it is in, i.e. the neutron generator's channel, and it's height via z).
 **/
void initialiseNeutron(struct neutron_struct * neutron, struct channel_struct * channel, double z) {
  neutron->active=true;
  generateNeutronVelocityAndEnergy(neutron);
  neutron->pos_x=channel->x_centre;
  neutron->pos_y=channel->y_centre;
  neutron->pos_z=z;
}

/**
 * Will calculate the amount of energy in MeV that a number of fissions
 * has produced
 **/
double getMeVFromFissions(unsigned long int num_fissions) {
  return num_fissions*MeV_PER_FISSION;
}

/**
 * Converts Mev to Joules
 **/
double getJoulesFromMeV(double MeV) {
  return JOULES_PER_MEV*MeV;
}

/**
 * Gets the number of atoms in a gram of a specific chemical
 **/
double getAtomsPerGram(enum chemical_type_enum chemical) {
  int restMass=getAtomicRestMass(chemical);
  return AVOGADRO/restMass;
}

/**
 * Given the simulations index of a chemical will convert this to
 * the element and return it
 **/
enum chemical_type_enum getChemicalAtIndex(int index) {
  if (index == 0) return U235;
  if (index == 1) return U238;
  if (index == 2) return Pu239;
  if (index == 3) return U236;
  if (index == 4) return Ba141;
  if (index == 5) return Kr92;
  if (index == 6) return Xe140;
  if (index == 7) return Sr94;
  if (index == 8) return Xe134;
  if (index == 9) return Zr103;
  if (index == 10) return Pu240;
  if (index == 11) return H2O;
  if (index == 12) return D2O;
  if (index == 13) return C6;
  return UNKNOWN_CHEMICAL;
}

/**
 * Determines whether a neutron scatters from the atom of the moderator based on its chemical composition,
 * the size of the moderator and its weight.
 *
 * This is a simplification, as the Barn unit is metres squared
 * whereas we assume all the atoms of the moderator are in that 2D plane, however realistically the neutron
 * is moving through the moderator so will likely come into contact with them in the other dimension anyway
 * in a timestep, so that is likely reasonable and it keeps things simpler too
 **/
static bool determineNeutronModeratorScattering(int moderator_weight, enum moderator_type_enum moderator_type, int size_z) {
  enum chemical_type_enum moderator_chemical;
  if (moderator_type == WATER) moderator_chemical=H2O;
  if (moderator_type == DEUTERIUM) moderator_chemical=D2O;
  if (moderator_type == GRAPHITE) moderator_chemical=C6;
  double cs=getScatteringCrossSection(moderator_chemical);

  double moderator_barns=size_z/BARN; // Number of barns of moderator height
  double size_reduction=moderator_barns/cs; // How many times larger the moderator is than the scattering cs

  // Get the number of atoms in the moderator
  double num_atoms=getAtomsPerGram(moderator_chemical) * moderator_weight;
  // Scale the number of atoms down by the size reduction to then see if there is likely one in that square
  double num_atoms_reduced=num_atoms/size_reduction;
  // If this number of atoms reduced is too small for the random number generator
  // then increase so there is always at-least a small chance
  double rand_resolution=getRandomNumberGeneratorResolution();
  if (num_atoms_reduced < rand_resolution) num_atoms_reduced=rand_resolution;
  double random_num=generateDecimalRandomNumber();
  return num_atoms_reduced >=random_num;
}

/**
 * Gets the scattering cross section of a specific chemical, here we are using
 * cross section constants that assume 1eV of energy of the neutron, which is
 * a simplification but it keeps the maths simpler
 **/
static double getScatteringCrossSection(enum chemical_type_enum chemical) {
  if (chemical == H2O) {
    return getScatteringCrossSectionOfWater();
  } else if (chemical == D2O) {
    return getScatteringCrossSectionOfDeuterium();
  } else if (chemical == C6) {
    return getScatteringCrossSectionOfGraphite();
  }
  return 0.0;
}

/**
 * Get the scattering cross section of water (H2O)
 **/
static double getScatteringCrossSectionOfWater() {
  int rest_mass=getAtomicRestMass(H2O);
  double n=(AVOGADRO/rest_mass) * DENSITY_WATER;
  return (2*n*SCATTERING_CROSS_SECTION_HYDROGEN) + (n*SCATTERING_CROSS_SECTION_OXYGEN);
}

/**
 * Get the scattering cross section of duterium (D2O)
 **/
static double getScatteringCrossSectionOfDeuterium() {
  int rest_mass=getAtomicRestMass(D2O);
  double n=(AVOGADRO/rest_mass) * DENSITY_DEUTERIUM;
  return (2*n*SCATTERING_CROSS_SECTION_HYDROGEN) + (n*SCATTERING_CROSS_SECTION_OXYGEN);
}

/**
 * Get the scattering cross section of graphite (C6)
 **/
static double getScatteringCrossSectionOfGraphite() {
  int rest_mass=getAtomicRestMass(C6);
  double n=(AVOGADRO/rest_mass) * DENSITY_GRAPHITE;
  return n*SCATTERING_CROSS_SECTION_CARBON;
}

/**
 * For a given neutron this will generate a random amount of energy in
 * MeV and the velocity components (which add up to 100 across the
 * three dimensions) that the velocity will be shared over.
 **/
static void generateNeutronVelocityAndEnergy(struct neutron_struct * neutron) {
  int sum=0;
  int x,y,z;
  while (sum == 0) {
    x=rand() % 100;
    y=rand() % 100;
    z=rand() % 100;
    sum = x+y+z;
  }
  neutron->x=x*(100/sum);
  if (rand() % 2 == 1) neutron->x=-neutron->x;
  neutron->y=y*(100/sum);
  if (rand() % 2 == 1) neutron->y=-neutron->y;
  neutron->z=z*(100/sum);
  if (rand() % 2 == 1) neutron->z=-neutron->z;
  // Random energy between 0 and 10 MeV
  neutron->energy=rand() / ((double) (RAND_MAX / 10));
}

/**
 * Determines whether a neutron has been absorbed by an atom of nuclear fuel within
 * a specific fuel pellet of a fuel rod.
 *
 * This is a simplification, as the Barn unit is metres squared
 * whereas we assume all the atoms of the fuel pellet are in that 2D plane, however
 * realistically the neutron is moving through the fuel channel so will likely come
 * into contact with them in the other dimension anyway
 * in a timestep, so that is likely reasonable and it keeps things simpler too
 **/
static bool determineNeutronAbsorbtionByFuel(double deBroglieWavelength, struct channel_struct * reactor_channel,
    enum chemical_type_enum chemical, int fuel_pellet, int collision_prob_multiplyer) {
  // First get the cross section in Barns for the chemical
  double cs=getCrossSectionInBarns(deBroglieWavelength, chemical);

  double pellet_b=0.0016/BARN; // Determine the size of the fuel pellet in barns
  double size_reduction=pellet_b/cs; // How much smaller the cross section is than the pellet

  // Now grab out the number of fuel atoms in the fuel pellet and scale this down by the cross section ratio
  double num_atoms=getNumberAtomsInFuelPellet(reactor_channel, chemical, fuel_pellet);
  double num_atoms_reduced=num_atoms/size_reduction;
  double rand_resolution=getRandomNumberGeneratorResolution();
  if (num_atoms_reduced < rand_resolution) num_atoms_reduced=rand_resolution;
  num_atoms_reduced *= collision_prob_multiplyer;
  // Now determine randomly if there is an atom in the cross section which will result in absorption
  double random_num=generateDecimalRandomNumber();
  return num_atoms_reduced >=random_num;
}

/**
 * Generates a random number between 0 and 1
 **/
static double generateDecimalRandomNumber() {
  return ((double) rand()) / RAND_MAX;
}

/**
 * Gets the resolution of the random number generator, this is important
 * when generating a random number between 0 and 1, as we might be comparing
 * against something that's so small it will never be chosen
 **/
static double getRandomNumberGeneratorResolution() {
  return 1.0/RAND_MAX;
}

/**
 * Retrieves the number of atoms of a chemical in a specific pellet of the fuel
 * held in the given reactor channel
 **/
static double getNumberAtomsInFuelPellet(struct channel_struct * reactor_channel,
    enum chemical_type_enum chemical, int fuel_pellet) {
  return reactor_channel->contents.fuel_assembly.quantities[fuel_pellet][chemical];
}

/**
 * Based upon a collision with a chemical, will calculate the cross section in Barnes from the
 * deBroglie wavelength that has been provided. 1 barn is 1e-28 meters squared.
 **/
static double getCrossSectionInBarns(double deBroglieWavelength, enum chemical_type_enum chemical) {
  int rest_mass=getAtomicRestMass(chemical);
  if (rest_mass == 0) {
    fprintf(stderr, "Rest mass of chemical determined as zero during cross section calculation\n");
    exit(-1);
  }
  double R=R0*pow(rest_mass, 0.333);
  double crossSection=PI*pow((R+(deBroglieWavelength/(2*PI))), 2);
  return crossSection/BARN;
}

/**
 * Returns the atomic rest mass of a chemical
 **/
static int getAtomicRestMass(enum chemical_type_enum chemical) {
  if (chemical == U235) return 235;
  if (chemical == U238) return 238;
  if (chemical == Pu239) return 239;
  if (chemical == U236) return 236;
  if (chemical == Ba141) return 141;
  if (chemical == Kr92) return 92;
  if (chemical == Xe140) return 140;
  if (chemical == Sr94) return 94;
  if (chemical == Xe134) return 134;
  if (chemical == Zr103) return 103;
  if (chemical == Pu240) return 240;
  if (chemical == H2O) return 18;
  if (chemical == D2O) return 19;
  if (chemical == C6) return 12;
  return 0;
}

/**
 * From the energy (in MeV) and atomic rest mass will calculate the
 * deBroglie wavelength
 **/
static double calculateDeBroglieWavelength(double MeV, int rest_mass) {
  double mass=rest_mass*MASS_ONE_UNIT;
  double ke=MeV*1e6*NEWTON_METER_PER_eV;
  double ppc=sqrt(ke*ke-(-1)*2*ke*mass*C*C);
  double pc=ppc/NEWTON_METER_PER_eV;
  double ddw=HC/pc;
  return ddw*NS_AS_SEC;
}
