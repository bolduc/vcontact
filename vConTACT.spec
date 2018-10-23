/*
A KBase module: vConTACT
*/

module vConTACT {
    /*
        Insert your typespec information here.
    */

    typedef structure {
        string parameter_1;
    } InParams;

    funcdef run_vcontact(InParams params)
        returns () authentication required;
};
