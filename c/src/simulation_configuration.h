#ifndef CONFIGURATION_INCLUDE
#define CONFIGURATION_INCLUDE

#include <stdbool.h>
#include "simulation_support.h"

// Maximum number of chemicals that will be configured
#define NUM_CONFIG_CHEMICALS 14
// Maximum number of reactor core channels in a row
#define MAX_CHANNEL_COLUMNS 1000
// Maximum number of control rods
#define MAX_CONTROL_RODS 1000

// Types of channel as described in configuration file
enum config_channel_type_enum {CONFIG_EMPTY=0, CONFIG_FUEL_ASSEMBLY=1, CONFIG_CONTROL_ROD=2, CONFIG_NEUTRON_GENERATOR=4, CONFIG_MODERATOR=5};
// Different chemicals that can be in the configuration file
enum config_chemical_type_enum {U235_CONFIG=0, U238_CONFIG=1, Pu239_CONFIG=2, U236_CONFIG=3, Ba141_CONFIG=4, Kr92_CONFIG=5, Xe140_CONFIG=6,
                                Sr94_CONFIG=7, Xe134_CONFIG=8, Zr103_CONFIG=9, Pu240_CONFIG=10, H2O_CONFIG=11, D2O_CONFIG=12, C6_CONFIG=13,
                                UNKNOWN_CHEMICAL_CONFIG=14};
// Type of moderator material
enum config_moderator_type_enum {WATER_MOD_TYPE_CONFIG, DEUTERIUM_MOD_TYPE_CONFIG, GRAPHITE_MOD_TYPE_CONFIG, NONE_MOD_TYPE_CONFIG};

// Holds configuration of control rod
struct control_rod_config_struct {
  int x_channel, y_channel, percentage;
};

// Overall configuration of the simulation
struct simulation_configuration_struct {
  int num_timesteps, dt, size_x, size_y, size_z, display_progess_frequency, write_reactor_state_frequency;
  int channels_x, channels_y, moderator_weight, collision_prob_multiplyer;
  enum config_channel_type_enum **channel_layout_config;
  enum config_moderator_type_enum moderator_type;
  int num_channel_configs[MAX_CHANNEL_COLUMNS];
  int fuel_makeup_percentage[NUM_CONFIG_CHEMICALS];
  struct control_rod_config_struct * control_rod_configurations;
  int num_ctrl_rod_configurations;
  long int max_neutrons;
};

void parseConfiguration(char*, struct simulation_configuration_struct*);
int findControlRodConfiguration(struct simulation_configuration_struct*, int, int);

#endif
