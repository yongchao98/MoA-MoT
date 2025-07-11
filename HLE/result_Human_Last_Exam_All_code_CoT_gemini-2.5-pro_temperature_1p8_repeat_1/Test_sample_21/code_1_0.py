import math

def calculate_density():
    """
    Calculates the density of primes for which f(x) is irreducible mod p
    by calculating the proportion of 7-cycles in candidate Galois groups.
    """
    
    # Candidate groups after preliminary analysis.
    # G must be a transitive subgroup of A_7, but not C_7.
    candidate_groups = {
        "F_21": { "order": 21 },
        "PSL(2,F_7)": { "order": 168 },
        "A_7": { "order": math.factorial(7) // 2 }
    }
    
    # --- Case 1: G = F_21 (Frobenius group of order 21) ---
    # In F_21 = C_7 semi-direct C_3, the elements of order 7 are the non-identity
    # elements of the normal C_7 subgroup.
    g1_name = "F_21"
    g1_order = candidate_groups[g1_name]["order"]
    # Number of elements of order 7 in C_7 is phi(7).
    g1_7_cycles = math.phi(7)
    g1_density_num = g1_7_cycles
    g1_density_den = g1_order
    
    # --- Case 2: G = PSL(2, F_7) (Simple group of order 168) ---
    # The number of Sylow 7-subgroups n_7 must be 1 mod 7 and divide |G|/7 = 24.
    # As PSL(2,F_7) is simple, n_7 > 1, so n_7 = 8.
    g2_name = "PSL(2,F_7)"
    g2_order = candidate_groups[g2_name]["order"]
    num_sylow_7_subgroups = 8
    # Each Sylow 7-subgroup has phi(7) elements of order 7.
    g2_7_cycles = num_sylow_7_subgroups * math.phi(7)
    g2_density_num = g2_7_cycles
    g2_density_den = g2_order
    
    # --- Case 3: G = A_7 (Alternating group on 7 elements) ---
    # The 7-cycles are even permutations. The number of n-cycles in S_n is (n-1)!.
    g3_name = "A_7"
    g3_order = candidate_groups[g3_name]["order"]
    g3_7_cycles = math.factorial(7 - 1)
    g3_density_num = g3_7_cycles
    g3_density_den = g3_order

    print("Calculating the density of primes p for which f(x) is irreducible mod p.")
    print("This is the proportion of 7-cycles in the Galois group G.\n")
    print("Based on properties of the polynomial, G must be F_21, PSL(2,F_7), or A_7.")
    print("Let's compute the proportion for each possibility:\n")
    
    print(f"Case 1: G = {g1_name}")
    print(f"Number of 7-cycles = {g1_7_cycles}")
    print(f"Order of group = {g1_order}")
    print(f"Proportion = {g1_density_num} / {g1_density_den} = {g1_density_num/math.gcd(g1_density_num, g1_density_den)}/{g1_density_den/math.gcd(g1_density_num, g1_density_den)}")
    print("-" * 20)
    
    print(f"Case 2: G = {g2_name}")
    print(f"Number of 7-cycles = {g2_7_cycles}")
    print(f"Order of group = {g2_order}")
    print(f"Proportion = {g2_density_num} / {g2_density_den} = {g2_density_num/math.gcd(g2_density_num, g2_density_den)}/{g2_density_den/math.gcd(g2_density_num, g2_density_den)}")
    print("-" * 20)

    print(f"Case 3: G = {g3_name}")
    print(f"Number of 7-cycles = {g3_7_cycles}")
    print(f"Order of group = {g3_order}")
    print(f"Proportion = {g3_density_num} / {g3_density_den} = {g3_density_num/math.gcd(g3_density_num, g3_density_den)}/{g3_density_den/math.gcd(g3_density_num, g3_density_den)}")
    print("-" * 20)

    final_num = g1_density_num//math.gcd(g1_density_num, g1_density_den)
    final_den = g1_density_den//math.gcd(g1_density_num, g1_density_den)
    
    print(f"\nSince the proportion is the same in all cases, the density is {final_num}/{final_den}.")

calculate_density()