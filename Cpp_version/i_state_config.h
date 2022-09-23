#ifndef I_STATE_CONFIG_H_
#define I_STATE_CONFIG_H_
#include "sp_util.h"

using namespace std;

const int max_id_num = 35000;
const int max_r_id_num = 2*max_id_num;

struct CFG {
    int max_r_num;
    string r_types[16];
    int pop_num;
    int total_num;
    vector<int> id_array;
    vector<int> r_id_array;
    vector<double> cr_vals;
    map<int,int> rid2id;
    map<int,string> rid2rtype;
    map<int,int> id2pid;
    double srv;
    double dt;
    double r1;
    double r2;
    int selected_r_id;
};

class Individual {
    public:
        struct i_state_vars {
            int ID;
            double age;
            string group;
        } isv;
        Individual() {
            isv.ID = 0;
            isv.age = 0.0;
            isv.group = "j";
        }; // constructor
        void i_based_reaction( int p_id, int id, vector<vector<int>> group_vals_array, vector<double> parm_vals, struct CFG& cfg ) {
    
            double r_rate[cfg.max_r_num];
            if( isv.group == "j" ) {
                r_rate[0] = 0.0;
                r_rate[1] = parm_vals[1]; // d_J
            } else if( isv.group == "a" ) {
                r_rate[0] = parm_vals[0]*exp(-1.0*group_vals_array[0][1]/parm_vals[3]);
                r_rate[1] = parm_vals[2]; // d_A
            }
    
            for( int i=0; i<cfg.max_r_num; i++ ) {
                cfg.srv += r_rate[i];
                cfg.cr_vals[cfg.max_r_num*id+i] = cfg.srv;
            }
        }
        //~Individual(); // destructor
};

#endif /* I_STATE_CONFIG_H_ */