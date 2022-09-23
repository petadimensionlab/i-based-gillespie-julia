//// MIT License ////
//The MIT License (MIT)
//Copyright (c) <2013> <Shinji Nakaoka>
//Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "individual_based_Gillespie_algorithm.cpp"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
//The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software. THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#pragma warning(disable:4786)
#include <iostream>
#include <cmath>
#include <ctime>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <random> // C++0x //
#include "ibm.h"
#include "i_state_config.h"

uniform_real_distribution<double> runif(0.0,1.0);

const double min_step = 0.1;
const double plot_interval = 1.0;
const double t0 = 0.0;
const double Tmax = 120.0;

void birth_death_process( vector<vector<Individual>>& population_array, vector<vector<int>>& group_vals_array, vector<map<string,int>> group_names_array, struct CFG& cfg ) {

    string selected_r_type = cfg.rid2rtype[cfg.selected_r_id];
    int selected_id = cfg.rid2id[cfg.selected_r_id];
    int selected_pid = cfg.id2pid[selected_id];
    if( selected_r_type == "+" ) {
        //  update state change by birth //
        int pop_number = get_pop_number(selected_pid,group_vals_array);
        int new_id = pop_number;
        Individual new_indiv;
        new_indiv = Individual();
        //new_indiv.isv.ID = new_id;
        population_array[selected_pid][new_id] = new_indiv;
        update_number(group_vals_array,selected_pid,selected_r_type,"j",group_names_array);
        update_configuration(selected_pid,selected_r_type,group_vals_array,cfg);
    } else if( selected_r_type == "-" ) {
        string prev_group = population_array[selected_pid][selected_id].isv.group;
        update_number(group_vals_array,selected_pid,selected_r_type,prev_group,group_names_array);
        population_array[selected_pid][selected_id] = population_array[selected_pid][cfg.total_num-1];
        //population_array[selected_pid][selected_id].isv.ID = selected_id;
        update_configuration(selected_pid,selected_r_type,group_vals_array,cfg);
    }
}

void maturation_process( vector<vector<Individual>>& population_array, vector<vector<int>>& group_vals_array, vector<double> parm_vals, vector<map<string,int>> group_names_array, double dt, struct CFG& cfg ) {

    int pop_number = get_pop_number(0,group_vals_array);
    for( int i=0; i<pop_number; i++ ) {
        double i_age = population_array[0][i].isv.age;
        string i_group = population_array[0][i].isv.group;
        population_array[0][i].isv.age += dt;
        if( i_age>=parm_vals[4] &&  i_group == "j" ) {
            population_array[0][i].isv.group = "a";
            update_number(group_vals_array,0,"+","a",group_names_array);
            update_number(group_vals_array,0,"-","j",group_names_array);
        }

    }
}

void initialize_ibm( vector<vector<Individual>>& population_array, vector<vector<int>> group_vals_array, vector<double> parm_vals, struct CFG& cfg ) {

    int p_count, g_count;
    Individual new_indiv;

    p_count = 0;
    for( int pid=0; pid<cfg.pop_num; pid++ ) {
        int m = group_vals_array[pid].size();
        g_count = 0; // reset g_count
        for( int id=0; id<m; id++ ) {
            int n = group_vals_array[pid][id];
            for( int k=0; k<n; k++ ) {
                new_indiv = Individual();
                new_indiv.isv.ID = p_count;
                new_indiv.isv.group = "a";
                population_array[pid][k] = new_indiv;
                set_initial_configuration(pid,p_count,cfg);
                p_count += 1;
                g_count += 1;
            }
        }
    }
    cfg.total_num = p_count;
    
}

int main(void) {

    int rnd_seed = 1;
    //random_device generate_seed;
    //mt19937 eng(generate_seed());
    mt19937 eng(rnd_seed);

    struct CFG cfg;
    cfg.total_num = 0;
    cfg.max_r_num = 2;
    cfg.pop_num = 1;
    const string pm[2] = {"+","-"};
    for( int i=0; i<2; i++ ) {
        cfg.r_types[i] = pm[i];
    }
    for( int i=0; i<max_id_num; i++ ) {
        cfg.id_array.push_back(-1);
    }
    for( int i=0; i<max_r_id_num; i++ ) {
        cfg.r_id_array.push_back(-1);
        cfg.cr_vals.push_back(0.0);
    }
    
    int parm_num = 5;
    vector<double> parm_vals(parm_num);
    parm_vals[0] = 8.5; // beta
    parm_vals[1] = 0.0060455567; // d_J
    parm_vals[2] = 0.27; // d_A
    parm_vals[3] = 600.0; // c
    parm_vals[4] = 15.6; // tau
    map<string,int> parm_names;
    parm_names.insert( pair<string,int>("beta",0) );
    parm_names.insert( pair<string,int>("dJ",1) );
    parm_names.insert( pair<string,int>("dA",2) );
    parm_names.insert( pair<string,int>("c",3) );
    parm_names.insert( pair<string,int>("tau",4) );

    vector<vector<int>> group_vals_array(cfg.pop_num);
    group_vals_array[0].push_back(0);
    group_vals_array[0].push_back(5000);
    vector<map<string,int>> group_names_array(cfg.pop_num);
    group_names_array[0].insert( pair<string,int>("j",0) );
    group_names_array[0].insert( pair<string,int>("a",1) );

    vector<Individual> population(max_id_num);
    vector<vector<Individual>> population_array(cfg.pop_num);
    for( int i=0; i<cfg.pop_num; i++ ) {
        population_array[i] = population;
    }
    initialize_ibm(population_array,group_vals_array,parm_vals,cfg);
    
    double t = t0;
    double next_plot = plot_interval;

    while( t<Tmax ) {
        //safety_check(population_array,group_vals_array,cfg);
        while( next_plot<t+plot_interval ) {
            if( next_plot>Tmax ) {
                break;
            }
            cout << "time = " << t << ", J = " << group_vals_array[0][0] << ", A = " << group_vals_array[0][1] << endl;
            next_plot += plot_interval;
        }
        i_based_Gillespie_direct(population_array,group_vals_array,parm_vals,eng,cfg);
        birth_death_process(population_array,group_vals_array,group_names_array,cfg);

        if( cfg.dt>min_step ) {
            double dt_tmp = -cfg.dt;
            double t_tmp = t;
            while( dt_tmp<0 ) {
                maturation_process(population_array,group_vals_array,parm_vals,group_names_array,min_step,cfg);
                dt_tmp += min_step;
                t_tmp += min_step;
            }
            // The resting part //
            maturation_process(population_array,group_vals_array,parm_vals,group_names_array,-dt_tmp,cfg);
        } else {
            maturation_process(population_array,group_vals_array,parm_vals,group_names_array,cfg.dt,cfg);
        }
        t += cfg.dt;
    }

    return 0;

}
