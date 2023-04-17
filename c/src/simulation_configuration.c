#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <ctype.h>
#include "simulation_configuration.h"

// Maximum length of configuration line
#define MAX_LINE_LENGTH 256

static void setControlRodConfiguration(struct simulation_configuration_struct*, char*);
static void checkConfiguration(struct simulation_configuration_struct*);
static void initialiseSimulationConfiguration(struct simulation_configuration_struct*);
static enum config_chemical_type_enum getChemicalEnum(char*);
static char* getChemicalName(char*);
static enum config_channel_type_enum* getChannelRowConfig(char*,int*);
static int countNumOccurances(char*, char);
static int getIntValue(char*);
static long getLongValue(char*);
static int getIntPercentValue(char*);
static int getEntityNumber(char*);
static double getDoubleValue(char*);
static char* getStringValue(char*);
static bool hasValue(char*);

/**
* A simple configuration file reader, I don't think you will need to change this (but feel free if you want to!)
* It will parse the configuration file and set the appropriate configuration points that will then feed into the simulation
* setup. It is somewhat limited in its flexibility and you need to be somewhat careful about the configuration file format,
* but is fine for our purposes
**/
void parseConfiguration(char * filename, struct simulation_configuration_struct * simulation_configuration) {
  initialiseSimulationConfiguration(simulation_configuration);
  FILE * f = fopen(filename, "r");
  if (f == NULL) {
    printf("Error, can not open file %s for reading\n", filename);
    exit(-1);
  }
  char buffer[MAX_LINE_LENGTH];
  while((fgets(buffer, MAX_LINE_LENGTH, f))!= NULL) {
    // If the string ends with a newline then remove this to make parsing simpler
    if (isspace(buffer[strlen(buffer)-1])) buffer[strlen(buffer)-1]='\0';
    if (strlen(buffer) > 0) {
      if (buffer[0]=='#') continue; // This line is a comment so ignore
      if (hasValue(buffer)) {
        if (strstr(buffer, "NUM_TIMESTEPS") != NULL) simulation_configuration->num_timesteps=getIntValue(buffer);
        if (strstr(buffer, "DT") != NULL) simulation_configuration->dt=getIntValue(buffer);
        if (strstr(buffer, "DISPLAY_PROGRESS_FREQUENCY") != NULL) simulation_configuration->display_progess_frequency=getIntValue(buffer);
        if (strstr(buffer, "WRITE_REACTOR_STATE_FREQUENCY") != NULL) simulation_configuration->write_reactor_state_frequency=getIntValue(buffer);
        if (strstr(buffer, "SIZE_X") != NULL) simulation_configuration->size_x=getIntValue(buffer);
        if (strstr(buffer, "SIZE_Y") != NULL) simulation_configuration->size_y=getIntValue(buffer);
        if (strstr(buffer, "SIZE_Z") != NULL) simulation_configuration->size_z=getIntValue(buffer);
        if (strstr(buffer, "MOD_WEIGHT") != NULL) simulation_configuration->moderator_weight=getIntValue(buffer);
        if (strstr(buffer, "COLLISION_PROB_MULTIPLYER") != NULL) simulation_configuration->collision_prob_multiplyer=getIntValue(buffer);
        if (strstr(buffer, "MAX_NEUTRONS") != NULL) simulation_configuration->max_neutrons=getLongValue(buffer);

        if (strstr(buffer, "CHANNELROW_") != NULL) {
          int channelRowNumber=getEntityNumber(buffer);
          if (channelRowNumber >= 0) {
            simulation_configuration->channel_layout_config[channelRowNumber]=getChannelRowConfig(buffer, &(simulation_configuration->num_channel_configs[channelRowNumber]));
          } else {
            fprintf(stderr, "Ignoring body configuration line '%s' as this is malformed and can not extract body number\n", buffer);
            exit(-1);
          }
        }

        if (strstr(buffer, "CONTROLROD_") != NULL) {
          setControlRodConfiguration(simulation_configuration, buffer);
        }

        if (strstr(buffer, "FUEL_") != NULL) {
          char * chemical=getChemicalName(buffer);
          if (chemical != NULL) {
            enum config_chemical_type_enum chemical_type=getChemicalEnum(chemical);
            if (chemical_type != UNKNOWN_CHEMICAL_CONFIG) {
              simulation_configuration->fuel_makeup_percentage[chemical_type]=getIntPercentValue(buffer);
            } else {
              fprintf(stderr, "Chemical not recognised for '%s'\n", chemical);
              exit(-1);
            }
            free(chemical);
          } else {
            fprintf(stderr, "Unknown chemical name in '%s'\n", buffer);
            exit(-1);
          }
        }

        if (strstr(buffer, "MODERATOR") != NULL) {
          char * moderator_type=getStringValue(buffer);
          if (strcmp(moderator_type, "WATER")==0) {
            simulation_configuration->moderator_type=WATER_MOD_TYPE_CONFIG;
          } else if (strcmp(moderator_type, "DEUTERIUM")==0) {
            simulation_configuration->moderator_type=DEUTERIUM_MOD_TYPE_CONFIG;
          } else if (strcmp(moderator_type, "GRAPHITE")==0) {
            simulation_configuration->moderator_type=GRAPHITE_MOD_TYPE_CONFIG;
          } else {
            simulation_configuration->moderator_type=NONE_MOD_TYPE_CONFIG;
          }
        }
      } else {
        fprintf(stderr, "Ignoring configuration line '%s' as this is malformed\n", buffer);
        exit(-1);
      }
    }
  }
  simulation_configuration->channels_x=(simulation_configuration->size_x * 100) / 20;
  simulation_configuration->channels_y=(simulation_configuration->size_y * 100) / 20;
  checkConfiguration(simulation_configuration);
  fclose(f);
}

/**
 * Finds the control rod configuration index for a specific channel in x and y
 **/
int findControlRodConfiguration(struct simulation_configuration_struct * simulation_configuration, int channel_x, int channel_y) {
  for (int i=0;i<simulation_configuration->num_ctrl_rod_configurations;i++) {
    if (simulation_configuration->control_rod_configurations[i].x_channel == channel_x &&
        simulation_configuration->control_rod_configurations[i].y_channel == channel_y) {
      return i;
    }
  }
  return -1;
}

/**
 * Sets the control rod at channel x and y as per _x_y to be inserted by a specified percentage
 **/
static void setControlRodConfiguration(struct simulation_configuration_struct * simulation_configuration, char * sourceString) {
  char * underScoreLocation=strchr(sourceString, '_');
  if (underScoreLocation != NULL) {
    char * secondUnderScoreLocation=strchr(underScoreLocation+1, '_');
    if (secondUnderScoreLocation != NULL) {
      char * equalsLocation=strchr(secondUnderScoreLocation+1, '=');
      int size_diff=secondUnderScoreLocation-underScoreLocation;
      int size_diff2=equalsLocation-secondUnderScoreLocation;
      char x_channel[size_diff], y_channel[size_diff2];
      strncpy(x_channel, &underScoreLocation[1], size_diff-1);
      x_channel[size_diff-1]='\0';
      strncpy(y_channel, &secondUnderScoreLocation[1], size_diff2-1);
      y_channel[size_diff2-1]='\0';
      simulation_configuration->control_rod_configurations[simulation_configuration->num_ctrl_rod_configurations].x_channel=atoi(x_channel);
      simulation_configuration->control_rod_configurations[simulation_configuration->num_ctrl_rod_configurations].y_channel=atoi(y_channel);
      simulation_configuration->control_rod_configurations[simulation_configuration->num_ctrl_rod_configurations].percentage=getIntPercentValue(sourceString);
      simulation_configuration->num_ctrl_rod_configurations++;
    } else {
      fprintf(stderr, "Control rod percentage configuration malformed '%s'\n", sourceString);
      exit(-1);
    }
  } else {
    fprintf(stderr, "Control rod percentage configuration malformed '%s'\n", sourceString);
    exit(-1);
  }
}

/**
 * Once configuration has been read then check it for validity before returning back to the main code, this
 * helps ensure that errornous values, or combinations of them, are not used in simulation
 **/
static void checkConfiguration(struct simulation_configuration_struct * simulation_configuration) {
  for (int i=0;i<NUM_CONFIG_CHEMICALS;i++) {
    // Check percentage of each chemical in fuel is a valid percentage value
    if (simulation_configuration->fuel_makeup_percentage[i] > 100 || simulation_configuration->fuel_makeup_percentage[i] < 0) {
      fprintf(stderr, "Reactor fuel provided at percentage '%d' but this is an incorrect percentage out of 100\n", simulation_configuration->fuel_makeup_percentage[i]);
      exit(-1);
    }
  }

  for (int i=0;i<simulation_configuration->num_ctrl_rod_configurations;i++) {
    // Check control rod configuration matches a control rod channel and percentage is valid
    int x_channel=simulation_configuration->control_rod_configurations[i].x_channel;
    int y_channel=simulation_configuration->control_rod_configurations[i].y_channel;
    if (simulation_configuration->control_rod_configurations[i].percentage > 100 ||
        simulation_configuration->control_rod_configurations[i].percentage < 0) {
      fprintf(stderr, "Control rod at x=%d y=%d specified at percentage '%d' but this is an incorrect percentage out of 100\n",
        x_channel, y_channel, simulation_configuration->control_rod_configurations[i].percentage);
      exit(-1);
    }
    if (x_channel == -1 || y_channel == -1) {
      fprintf(stderr, "Active control rod configuration has x=%d y=%d where -1 is uninitialised\n", x_channel, y_channel);
      exit(-1);
    }
    if (simulation_configuration->channel_layout_config[x_channel] != NULL) {
      if (simulation_configuration->num_channel_configs[x_channel] > y_channel) {
        if (simulation_configuration->channel_layout_config[x_channel][y_channel] != CONFIG_CONTROL_ROD) {
          fprintf(stderr, "Control rod configuration provided for channel x=%d y=%d but this type is not a control rod\n", x_channel, y_channel);
          exit(-1);
        }
      } else {
        fprintf(stderr, "Control rod configuration provided for channel x=%d y=%d but this is an empty channel\n", x_channel, y_channel);
        exit(-1);
      }
    } else {
      fprintf(stderr, "Control rod configuration provided for channel x=%d y=%d but this is an empty channel\n", x_channel, y_channel);
      exit(-1);
    }
  }

  for (int i=0;i<MAX_CHANNEL_COLUMNS;i++) {
    // Check channels that are configured are sensible
    if (simulation_configuration->channel_layout_config[i] != NULL) {
      if (i >= simulation_configuration->channels_x) {
        fprintf(stderr, "Reactor channel row %d configured, but there are only %d rows\n", i, simulation_configuration->channels_x);
        exit(-1);
      }
      if (simulation_configuration->num_channel_configs[i] > simulation_configuration->channels_y) {
        fprintf(stderr, "There are %d reactor channel columns, but row %d has %d\n", simulation_configuration->channels_y,
          i, simulation_configuration->num_channel_configs[i]);
        exit(-1);
      }
    }
  }
}

/**
* Initialises the configuration data structure for this simulation, adds in some defaults if they are not specified and marks
* each channel as inactive
**/
static void initialiseSimulationConfiguration(struct simulation_configuration_struct * simulation_configuration) {
  // Default values
  simulation_configuration->dt=10;
  simulation_configuration->num_timesteps=1000;
  simulation_configuration->num_ctrl_rod_configurations=0;
  simulation_configuration->write_reactor_state_frequency=0;
  simulation_configuration->display_progess_frequency=0;
  simulation_configuration->channel_layout_config=(enum config_channel_type_enum **) malloc(sizeof(enum config_channel_type_enum *) * MAX_CHANNEL_COLUMNS);
  simulation_configuration->control_rod_configurations=(struct control_rod_config_struct *) malloc(sizeof(struct control_rod_config_struct) * MAX_CONTROL_RODS);
  for (int i=0;i<MAX_CHANNEL_COLUMNS;i++) {
    simulation_configuration->channel_layout_config[i]=NULL;
    simulation_configuration->num_channel_configs[i]=0;
  }
  for (int i=0;i<MAX_CONTROL_RODS;i++) {
    simulation_configuration->control_rod_configurations[i].x_channel=-1;
    simulation_configuration->control_rod_configurations[i].y_channel=-1;
  }
  for (int i=0;i<NUM_CONFIG_CHEMICALS;i++) {
    simulation_configuration->fuel_makeup_percentage[i]=0;
  }
}

/**
 * Mapping between the name of chemicals provided in the configuration file and their configuration enumeration
 **/
static enum config_chemical_type_enum getChemicalEnum(char * name) {
  if (strcmp(name, "U235")==0) return U235_CONFIG;
  if (strcmp(name, "U238")==0) return U238_CONFIG;
  if (strcmp(name, "Pu239")==0) return Pu239_CONFIG;
  if (strcmp(name, "U236")==0) return U236_CONFIG;
  if (strcmp(name, "Ba141")==0) return Ba141_CONFIG;
  if (strcmp(name, "Kr92")==0) return Kr92_CONFIG;
  if (strcmp(name, "Xe140")==0) return Xe140_CONFIG;
  if (strcmp(name, "Sr94")==0) return Sr94_CONFIG;
  if (strcmp(name, "Xe134")==0) return Xe134_CONFIG;
  if (strcmp(name, "Zr103")==0) return Zr103_CONFIG;
  if (strcmp(name, "Pu240")==0) return Pu240_CONFIG;
  return UNKNOWN_CHEMICAL_CONFIG;
}

/**
 * Retrieves the chemical name that has been added after the underscore,
 * for instance FUEL_U235 will return U235
 **/
static char* getChemicalName(char * sourceString) {
  char * underScoreLocation=strchr(sourceString, '_');
  if (underScoreLocation != NULL) {
    char * secondUnderScoreLocation=strchr(underScoreLocation+1, '=');
    if (secondUnderScoreLocation != NULL) {
      int size_diff=secondUnderScoreLocation-underScoreLocation;
      char * chemical=(char*) malloc(size_diff);
      strncpy(chemical, &underScoreLocation[1], size_diff-1);
      chemical[size_diff-1]='\0';
      return chemical;
    }
  }
  return NULL;
}

/**
 * Based upon a source string will figure out the reactor core's channels configurations
 * and return this, along with the number of options specified. For instance
 * "FUEL,CONTROL,MODERATOR,CONTROL,FUEL" would return enumerations matching those
 * strings and five as the number of options
 **/
static enum config_channel_type_enum* getChannelRowConfig(char * sourceString, int * number_options_specified) {
  char cmp_string[50];
  char * equalsLocation=strchr(sourceString, '=');
  if (equalsLocation != NULL) {
    *number_options_specified=countNumOccurances(equalsLocation+1, ',')+1;
    enum config_channel_type_enum* row_layout_config=(enum config_channel_type_enum*) malloc(sizeof(enum config_channel_type_enum) * *number_options_specified);
    char * searchString=equalsLocation+1;
    for (int i=0;i<*number_options_specified;i++) {
      char * comma_Location;
      if (i < *number_options_specified-1) {
        comma_Location=strchr(searchString, ',');
      } else {
        comma_Location=searchString+strlen(searchString);
      }
      int option_len=comma_Location-searchString;
      strncpy(cmp_string, searchString, option_len);
      cmp_string[option_len]='\0';
      if (strcmp(cmp_string, "FUEL")==0) {
        row_layout_config[i]=CONFIG_FUEL_ASSEMBLY;
      } else if (strcmp(cmp_string, "CONTROL")==0) {
        row_layout_config[i]=CONFIG_CONTROL_ROD;
      } else if (strcmp(cmp_string, "NEUTRON")==0) {
        row_layout_config[i]=CONFIG_NEUTRON_GENERATOR;
      } else if (strcmp(cmp_string, "MODERATOR")==0) {
        row_layout_config[i]=CONFIG_MODERATOR;
      } else if (strcmp(cmp_string, "EMPTY")==0) {
        row_layout_config[i]=CONFIG_EMPTY;
      }
      searchString=comma_Location+1;
    }
    return row_layout_config;
  }
  return NULL;
}

/**
 * Helper to count the number of occurances of the character 'c' in
 * a string that has been provided
 **/
static int countNumOccurances(char * sourceString, char c) {
  int str_len=strlen(sourceString);
  int num_char=0;
  for (int i=0;i<str_len;i++) {
    if (sourceString[i] == c) num_char++;
  }
  return num_char;
}

/**
* From a string will extract the integer value after the equals and
* return this, or -1 if none is found
**/
static int getIntValue(char * sourceString) {
  char * equalsLocation=strchr(sourceString, '=');
  if (equalsLocation != NULL) {
    return atoi(&equalsLocation[1]);
  }
  return -1;
}

/**
* From a string will extract a long integer value after the
* equals and return this, or -1 if none is found
**/
static long getLongValue(char * sourceString) {
  char * equalsLocation=strchr(sourceString, '=');
  if (equalsLocation != NULL) {
    return atol(&equalsLocation[1]);
  }
  return -1;
}

/**
* From a string will extract the percentage integer value after the
* equals and return this, or -1 if none is found
**/
static int getIntPercentValue(char * sourceString) {
  char * equalsLocation=strchr(sourceString, '=');
  if (equalsLocation != NULL) {
    char * pcLocation=strchr(equalsLocation+1, '%');
    if (pcLocation == NULL) {
      return atoi(&equalsLocation[1]);
    } else {
      int size_diff=pcLocation-equalsLocation;
      char int_val[size_diff];
      strncpy(int_val, &equalsLocation[1], size_diff-1);
      int_val[size_diff-1]='\0';
      return atoi(int_val);
    }
  }
  return -1;
}

/**
* A helper function to parse a string with an underscore in it, this will extract the number after the underscore
* as we use this in the configuration file in a number of places
**/
static int getEntityNumber(char * sourceString) {
  char * underScoreLocation=strchr(sourceString, '_');
  if (underScoreLocation != NULL) {
    char * secondUnderScoreLocation=strchr(underScoreLocation+1, '=');
    if (secondUnderScoreLocation != NULL) {
      int size_diff=secondUnderScoreLocation-underScoreLocation;
      char int_key[size_diff];
      strncpy(int_key, &underScoreLocation[1], size_diff-1);
      int_key[size_diff-1]='\0';
      return atoi(int_key);
    }
  }
  return -1;
}

/**
* From a string will extract the double value after the equals and return
* this, or -1 if none is found
**/
static double getDoubleValue(char * sourceString) {
  char * equalsLocation=strchr(sourceString, '=');
  if (equalsLocation != NULL) {
    return atof(&equalsLocation[1]);
  }
  return -1;
}

/**
* From a string will extract the string after the equals and return this,
* or -1 if none is found
**/
static char * getStringValue(char * sourceString) {
  char * equalsLocation=strchr(sourceString, '=');
  if (equalsLocation != NULL) {
    return equalsLocation+1;
  }
  return NULL;
}

/**
* Determines if a string has an equals in it or not
* (e.g. is there a value specified at this line?)
**/
static bool hasValue(char * sourceString) {
  return strchr(sourceString, '=') != NULL;
}
