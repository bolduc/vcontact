/*
A KBase module: vConTACT
*/

module vConTACT {
    /*
        Insert your typespec information here.

        **ERROR** - KIDL-spec validation failed:
                Error at line 32, column 9: Can not find type: bool
    */

    typedef string obj_ref;

    typedef structure {
        string workspace_name;
        obj_ref genome;
        string db;
        string pcs_mode;
        string vcs_mode;
        float blast_evalue;
        float pc_inflation;
        float pc_max_overlap;
        float pc_penalty;
        float pc_haircut;
        float vc_inflation;
        float vc_density;
        int vc_min_size;
        float vc_max_overlap;
        int vc_penalty;
        float vc_haircut;
        string merge_method;
        string similarity;
        string seed_method;
        string optimize;
        float min_significance;
        float max_significance;
        string permissive;
        float module_inflation;
        float mod_significance;
        int module_min_shared;
        float link_significance;
        float link_proportion;

    } InParams;

    funcdef run_vcontact(InParams params)
        returns () authentication required;
};
