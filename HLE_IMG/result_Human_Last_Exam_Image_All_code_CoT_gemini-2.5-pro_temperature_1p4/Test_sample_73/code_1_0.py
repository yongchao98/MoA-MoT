def get_stereochemical_assignments():
    """
    This function determines and prints the stereochemical assignments for the four stereocenters
    in the reaction scheme, from left to right.
    """
    # The four stereochemical assignments are determined as follows:
    # 1. The first reactant (acyl chloride): (S)
    # 2. The second reactant (alcohol): (S)
    # 3. The first stereocenter of the product (from the acyl chloride): (S)
    # 4. The second stereocenter of the product (from the alcohol): (S)
    assignments = ["S", "S", "S", "S"]
    print(", ".join(assignments))

get_stereochemical_assignments()