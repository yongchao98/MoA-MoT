def assign_stereochemistry():
    """
    Determines and prints the stereochemical assignments for the four stereocenters
    in the provided reaction scheme.
    """
    
    # Based on CIP rules and 3D structure analysis:
    # 1. The first stereocenter (acyl chloride reactant).
    assignment_1 = "(R)"
    
    # 2. The second stereocenter (alcohol reactant).
    assignment_2 = "(R)"
    
    # The reaction is an esterification which does not affect the stereocenters.
    # 3. The third stereocenter (in the product, from the alcohol).
    assignment_3 = "(R)" # Same as the second
    
    # 4. The fourth stereocenter (in the product, from the acyl chloride).
    assignment_4 = "(R)" # Same as the first
    
    print("Stereochemical assignments moving from left to right:")
    print(f"1. Acyl Chloride Reactant Stereocenter: {assignment_1}")
    print(f"2. Alcohol Reactant Stereocenter: {assignment_2}")
    print(f"3. Product Stereocenter (alcohol-derived): {assignment_3}")
    print(f"4. Product Stereocenter (acid-derived): {assignment_4}")

assign_stereochemistry()