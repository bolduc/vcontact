/*
A KBase module: vConTACT
*/

module vConTACT {
    /*
        Insert your typespec information here.
    */

    typedef string obj_ref;

    typedef structure {
        string workspace_name;
        obj_ref genome;
    } InParams;

    funcdef run_vcontact(InParams params)
        returns () authentication required;
};
