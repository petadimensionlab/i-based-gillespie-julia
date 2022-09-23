#include <iostream>
#include <cmath>
#include <ctime>
#include <string>
#include <vector>
#include <map>
#include <random>
#include <algorithm>
#include "i_state_config.h"

using namespace std;
extern uniform_real_distribution<double> runif;
extern const double min_step;

void generate_random_numbers( mt19937& eng, struct CFG& cfg );
void update_number( vector<vector<int>>& group_vals_array, int p_id, string r_type, string group_name, vector<map<string,int>> group_names_array );
int get_total_number( vector<vector<int>> group_vals_array );
int get_pop_number( int pid, vector<vector<int>> group_vals_array );
int get_cmr_pop_number( int pid, vector<vector<int>> group_vals_array );
void safety_check( vector<vector<Individual>> population_array, vector<vector<int>> group_vals_array, struct CFG& cfg );

void pid_shift2right( int pid, vector<vector<int>>& group_vals_array, struct CFG& cfg );
void pid_shift2left( int pid, vector<vector<int>>& group_vals_array, struct CFG& cfg );
void update_configuration( int selected_pid, string selected_r_type, vector<vector<int>>& group_vals_array, struct CFG& cfg );
void set_initial_configuration( int pid, int p_count, struct CFG& cfg );

void generate_random_numbers( mt19937& eng, struct CFG& cfg );
void i_based_Gillespie_direct( vector<vector<Individual>> population_array, vector<vector<int>> group_vals_array, vector<double> parm_vals, mt19937& eng, struct CFG& cfg );

