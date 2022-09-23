#include "ibm.h"

void update_number( vector<vector<int>>& group_vals_array, int p_id, string r_type, string group_name, vector<map<string,int>> group_names_array ) {
    int group_id;
    group_id = group_names_array[p_id][group_name];
    if( r_type == "+" ) {
        group_vals_array[p_id][group_id] += 1;
    } else if( r_type == "-" ) {
        group_vals_array[p_id][group_id] += -1;
    }
}

int get_total_number( vector<vector<int>> group_vals_array ) {
    int total_number = 0;
    int m = group_vals_array.size();
    for( int i=0; i<m; i++ ) {
        int n = group_vals_array[i].size();
        for( int j=0; j<n; j++ ) {
            total_number += group_vals_array[i][j];
        }
    }
    return total_number;
}

int get_pop_number( int pid, vector<vector<int>> group_vals_array ) {
    int pop_number = 0;
    int n = group_vals_array[pid].size();
    for( int j=0; j<n; j++ ) {
        pop_number += group_vals_array[pid][j];
    }
    return pop_number;
}

int get_cmr_pop_number( int pid, vector<vector<int>> group_vals_array ) {
    
    int pop_number = 0;
    for( int p=0; p<=pid; p++ ) {
        int n = group_vals_array[p].size();
        for( int j=0; j<n; j++ ) {
            pop_number += group_vals_array[p][j];
        }
    }
    return pop_number;
}
    
void safety_check( vector<vector<Individual>> population_array, vector<vector<int>> group_vals_array, struct CFG& cfg ) {
    int total_number = get_total_number(group_vals_array);
    if( total_number==0 ) {
        cout << "No individual is found. Stop computation." << endl;
        exit(0);
    }
    for( int pid; pid<cfg.pop_num; pid++ ) {
        int pop_number = get_pop_number(pid,group_vals_array);
        if( pop_number==0 ) {
            cout << "No individual is found in population " << pid << ". Stop computation." << endl;
            exit(0);
        }
    }
    
    if( cfg.total_num>max_id_num ) {
        cfg.id_array.resize(2*max_id_num);
        cfg.r_id_array.resize(2*max_r_id_num);
        cfg.cr_vals.resize(2*max_r_id_num);
        for( int i=0; i<cfg.pop_num; i++ ) {
            population_array[i].resize(2*max_id_num);
        }
    }
    if( std::isinf(cfg.dt)==1 ) {
        cout << "'dt' is infinity. Stop computation." << endl;
        exit(1);
    }
}

void pid_shift2right( int pid, vector<vector<int>>& group_vals_array, struct CFG& cfg ) {
    
    int cmr_pop_number;
    if( pid == cfg.pop_num-1 ) {
        // not necessary to shift //
        cfg.id2pid[cfg.total_num] = pid; // add a new item at the end | old: cfg.total_num+1
    } else {
        // shift 
        for( int p=pid; p<cfg.pop_num-1; p++ ) {
            cmr_pop_number = get_cmr_pop_number(p,group_vals_array);
            cfg.id2pid[cmr_pop_number+1] = p+1;
        }
        cfg.id2pid[cfg.total_num] = cfg.pop_num-1; // add a new item at the end | old: cfg.total_num+1
    }
    
}

void pid_shift2left( int pid, vector<vector<int>>& group_vals_array, struct CFG& cfg ) {
    
    int cmr_pop_number;
    if( pid == cfg.pop_num-1 ) {
        // not necessary to shift //
        cfg.id2pid[cfg.total_num] = -1;
    } else {
        // shift 
        for( int p=pid; p<cfg.pop_num-1; p++ ) {
            cmr_pop_number = get_cmr_pop_number(p,group_vals_array);
            cfg.id2pid[cmr_pop_number-1] = p+1;
        }
    }  
}

void update_configuration( int selected_pid, string selected_r_type, vector<vector<int>>& group_vals_array, struct CFG& cfg ) {

    if( selected_r_type == "+" ) {
        cfg.id_array[cfg.total_num] = cfg.total_num;
        pid_shift2right(selected_pid,group_vals_array,cfg);
        for( int i=0; i<cfg.max_r_num; i++ ) {
            cfg.r_id_array[cfg.max_r_num*cfg.total_num+i] = cfg.max_r_num*cfg.id_array[cfg.total_num]+i;
            cfg.rid2id[cfg.total_num*cfg.max_r_num+i] = cfg.total_num;
            cfg.rid2rtype[cfg.total_num*cfg.max_r_num+i] = cfg.r_types[i];
        }
        cfg.total_num += 1;
    } else if ( selected_r_type == "-" ) {
        pid_shift2left(selected_pid,group_vals_array,cfg);
        cfg.total_num += -1;
    }

}

void set_initial_configuration( int pid, int p_count, struct CFG& cfg ) {

    cfg.id_array[p_count] = p_count;
    cfg.id2pid[p_count] = pid;
    for( int i=0; i<cfg.max_r_num; i++ ) {
        cfg.r_id_array[cfg.max_r_num*p_count+i] = cfg.max_r_num*cfg.id_array[p_count]+i;
        cfg.rid2id[p_count*cfg.max_r_num+i] = p_count;
        cfg.rid2rtype[p_count*cfg.max_r_num+i] = cfg.r_types[i];
    }
}

void generate_random_numbers( mt19937& eng, struct CFG& cfg ) {

    cfg.r1 = runif(eng);
    cfg.r2 = runif(eng);
}

void i_based_Gillespie_direct( vector<vector<Individual>> population_array, vector<vector<int>> group_vals_array, vector<double> parm_vals, mt19937& eng, struct CFG& cfg ) {
    
    generate_random_numbers(eng,cfg);
    cfg.srv = 0.0;
    int cmr_count = 0;

    for( int pid=0; pid<cfg.pop_num; pid++ ) {
        int n = group_vals_array[pid].size();
        int group_count = 0;
        for( int j=0; j<n; j++ ) {
            group_count += group_vals_array[pid][j];
        }
        for( int count=0; count<group_count; count++ ) {
            population_array[pid][count].i_based_reaction(pid,cmr_count,group_vals_array,parm_vals,cfg);
            cmr_count += 1;
        }
    }
    
    if( cfg.srv==0.0 ) {
        cfg.dt = min_step;
        cfg.selected_r_id = -1;
    } else {
        double val;
        double thr = cfg.r2*cfg.srv;
        
        cfg.selected_r_id = 0;
        int i = 1;
        val = cfg.cr_vals[0];
        if( val>=thr ) {
            cfg.selected_r_id = cfg.r_id_array[0];
        } else {
            while( val<thr ) {
                val = cfg.cr_vals[i];
                i += 1;
            }
            cfg.selected_r_id = cfg.r_id_array[i-1];
        }
        cfg.dt = -log(cfg.r1)/cfg.srv;
    }
}


