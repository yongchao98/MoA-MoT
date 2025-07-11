def calculate_ni_ce_ratio():
    """
    Calculates and explains the ideal Ni/Ce atomic ratio based on literature findings.
    
    The ideal ratio for Ni-Ceria catalysts is often expressed in terms of atomic percent (at%).
    A value of 15 at% Ni is frequently found to be highly effective. This function
    calculates the corresponding Ni/Ce atomic ratio.
    """
    
    # A representative optimal value for Ni concentration in atomic percent (at%).
    ni_atomic_percent = 15
    
    # The remainder is the cerium (Ce) concentration.
    ce_atomic_percent = 100 - ni_atomic_percent
    
    # Calculate the atomic ratio of Ni to Ce.
    ni_ce_ratio = ni_atomic_percent / ce_atomic_percent
    
    print("The ideal Ni/Ce ratio is not a single value but falls within a range.")
    print(f"Based on research, a representative optimal composition is around {ni_atomic_percent} atomic percent Ni.")
    print("This means for every 100 total metal atoms (Ni + Ce), there are 15 Ni atoms and 85 Ce atoms.")
    print("\nThe corresponding Ni/Ce atomic ratio is calculated as follows:")
    print(f"Equation: Ni atoms / Ce atoms = {ni_atomic_percent} / {ce_atomic_percent}")
    print(f"Resulting Ratio: {ni_ce_ratio:.3f}")

calculate_ni_ce_ratio()