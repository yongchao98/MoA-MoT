def get_stereochemical_assignments():
    """
    This function returns the stereochemical assignments for the four stereocenters
    in the reaction scheme from left to right.
    """
    assignments = ["S", "S", "R", "S"]
    print("The stereochemical assignments for the four stereocenters from left to right are:")
    print(f"1st (Acyl Chloride): {assignments[0]}")
    print(f"2nd (Alcohol): {assignments[1]}")
    print(f"3rd (Product center 1): {assignments[2]}")
    print(f"4th (Product center 2): {assignments[3]}")
    print(f"\nFinal sequence: {', '.join(assignments)}")

get_stereochemical_assignments()