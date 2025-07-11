import math

def format_group(free_rank, torsion_parts):
    """Formats an abelian group into the specified string format."""
    parts = []
    if free_rank > 0:
        if free_rank == 1:
            parts.append("Z")
        else:
            parts.append(f"Z^{free_rank}")
    
    for t in sorted(torsion_parts):
        parts.append(f"Z/{t}Z")
        
    if not parts:
        return "0"
    
    return "+".join(parts)

def combinations(n, k):
    """Calculates the binomial coefficient 'n choose k'."""
    if k < 0 or k > n:
        return 0
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

def calculate_cohomology_M7():
    """
    Calculates the cohomology groups of M(7) based on known and assumed homology.
    """
    k = 7
    
    # Step 1: Calculate Homology Groups H_i(M(7))
    # H_0
    h0_free = 1
    h0_torsion = []
    
    # H_1
    h1_free = combinations(k, 2)
    h1_torsion = [2]
    
    # H_2
    h2_free_num = k * (k - 1) * (k - 2) * (3 * k - 5)
    h2_free_den = 24
    h2_free = h2_free_num // h2_free_den
    h2_torsion = [2, 2]
    
    # H_i for i >= 3 are assumed to be 0
    h3_free = 0
    h3_torsion = []

    homology = [
        (h0_free, h0_torsion),
        (h1_free, h1_torsion),
        (h2_free, h2_torsion),
        (h3_free, h3_torsion)
    ]

    # Step 2: Apply Universal Coefficient Theorem to get Cohomology H^i(M(7))
    cohomology_groups = []
    
    # H^0
    # Hom(H_0, Z) = Hom(Z, Z) = Z
    # Ext(H_{-1}, Z) = 0
    h_up_0_free = h0_free 
    h_up_0_torsion = []
    cohomology_groups.append(format_group(h_up_0_free, h_up_0_torsion))
    
    # H^1
    # Hom(H_1, Z) = Hom(Z^21 + Z/2Z, Z) = Z^21
    # Ext(H_0, Z) = Ext(Z, Z) = 0
    h_up_1_free = h1_free
    h_up_1_torsion = [] # from Ext(H_0, Z)
    cohomology_groups.append(format_group(h_up_1_free, h_up_1_torsion))

    # H^2
    # Hom(H_2, Z) = Hom(Z^140 + (Z/2Z)^2, Z) = Z^140
    # Ext(H_1, Z) = Ext(Z^21 + Z/2Z, Z) = Z/2Z
    h_up_2_free = h2_free
    h_up_2_torsion = h1_torsion # from Ext(H_1, Z)
    cohomology_groups.append(format_group(h_up_2_free, h_up_2_torsion))
    
    # H^3
    # Hom(H_3, Z) = Hom(0, Z) = 0
    # Ext(H_2, Z) = Ext(Z^140 + (Z/2Z)^2, Z) = (Z/2Z)^2
    h_up_3_free = h3_free
    h_up_3_torsion = h2_torsion # from Ext(H_2, Z)
    cohomology_groups.append(format_group(h_up_3_free, h_up_3_torsion))
    
    # H^i for i >= 4 will be 0 under our assumption
    
    print(f"H^0(M(7)) = {cohomology_groups[0]}")
    print(f"H^1(M(7)) = {cohomology_groups[1]} (where the rank is {h1_free})")
    print(f"H^2(M(7)) = {cohomology_groups[2]} (where the rank is {h2_free})")
    print(f"H^3(M(7)) = {cohomology_groups[3]}")
    print(f"\nThe list of cohomology groups is:")
    print(f"[{', '.join(cohomology_groups)}]")

calculate_cohomology_M7()