#ifndef SUPPORT_INCLUDE
#define SUPPORT_INCLUDE

#include <stdbool.h>

// Total number of chemicals that the simulation can represent
#define NUM_CHEMICALS 14

// Type of reactor channel
enum channel_type_enum {MODERATOR=0, FUEL_ASSEMBLY=1, CONTROL_ROD=2, EMPTY_CHANNEL=3, NEUTRON_GENERATOR=4};
// Type of moderator material
enum moderator_type_enum {WATER, DEUTERIUM, GRAPHITE, NO_MODERATOR};
// Chemicals that can be represented in the simulation
enum chemical_type_enum {U235=0, U238=1, Pu239=2, U236=3, Ba141=4, Kr92=5, Xe140=6, Sr94=7, Xe134=8, Zr103=9, Pu240=10, H2O=11, D2O=12, C6=13, UNKNOWN_CHEMICAL=14};

// Holds the state of the moderator channel
struct moderator_struct {
  enum moderator_type_enum type;
};

// Holds the state of a fuel assembly in a channel, fuel assemblies
// are made up of fuel pellets (40mm by 40mm by 2mm) stacked on top of each other
struct fuel_assembly_struct {
  double (*quantities)[NUM_CHEMICALS];
  int num_pellets;
  unsigned long int num_fissions;
};

// State of a control rod in a channel
struct control_rod_struct {
  double lowered_to_level;
};

// State of a neutron generator in a channel
struct neutron_generator_struct {
  double weight;
};

// Represents an individual neutron
struct neutron_struct {
  short x,y,z;
  double pos_x, pos_y, pos_z;
  double energy;
  bool active;
};

// A single channel in the reactor
struct channel_struct {
  union {
    struct moderator_struct moderator;
    struct fuel_assembly_struct fuel_assembly;
    struct control_rod_struct control_rod;
    struct neutron_generator_struct neutron_generator;
  } contents;

  enum channel_type_enum type;
  double x_centre, y_centre;
};

unsigned long int getNumberNeutronsFromGenerator(double, int);
double MeVToVelocity(double, int);
bool determineAndHandleIfNeutronModeratorCollision(struct neutron_struct*, int, enum moderator_type_enum, int);
bool determineAndHandleIfNeutronFuelCollision(double, struct channel_struct*, int, int);
int fissionPu240(struct channel_struct*, int);
int fissionU236(struct channel_struct*, int);
void initialiseNeutron(struct neutron_struct*, struct channel_struct*, double);
double getMeVFromFissions(unsigned long int);
double getJoulesFromMeV(double);
double getAtomsPerGram(enum chemical_type_enum);
enum chemical_type_enum getChemicalAtIndex(int);

#endif
