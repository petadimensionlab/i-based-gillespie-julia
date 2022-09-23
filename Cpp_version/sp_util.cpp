#include "sp_util.h"

string int2string( int number ) {

    string s;
    stringstream out;
    out << number;
    s = out.str();
    return s;
}

int string2int( string s ) {

    int number;
    number = (int)atoi(s.c_str());
    return number;
}
