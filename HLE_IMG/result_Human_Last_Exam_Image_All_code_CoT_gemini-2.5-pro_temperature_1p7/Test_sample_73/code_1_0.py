def solve_stereochemistry():
    """
    This function stores and prints the stereochemical assignments for the reaction.
    """
    # Assignments for the four stereocenters from left to right in the scheme.
    # 1: Acyl chloride reactant
    # 2: Alcohol reactant
    # 3: Product center from acyl chloride
    # 4: Product center from alcohol
    assignments = ["S", "S", "R", "S"]

    print("The stereochemical assignments for the four stereocenters are:")
    print(f"Stereocenter 1 (acyl chloride): ({assignments[0]})")
    print(f"Stereocenter 2 (alcohol): ({assignments[1]})")
    print(f"Stereocenter 3 (product, left): ({assignments[2]})")
    print(f"Stereocenter 4 (product, right): ({assignments[3]})")
    print("\nThe sequence from left to right is:")
    print(f"({assignments[0]}, {assignments[1]}, {assignments[2]}, {assignments[3]})")

solve_stereochemistry()