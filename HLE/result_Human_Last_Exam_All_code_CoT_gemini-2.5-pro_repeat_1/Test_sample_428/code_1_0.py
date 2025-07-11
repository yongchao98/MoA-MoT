import sys

def solve():
    """
    Calculates the expected time until only one particle remains.
    The user should provide the four integer values N1, M1, N2, M2 as command-line arguments.
    """
    if len(sys.argv) != 5:
        print("Usage: python your_script_name.py N1 M1 N2 M2")
        print("Example: python your_script_name.py 10 20 15 25")
        return

    try:
        N1 = int(sys.argv[1])
        M1 = int(sys.argv[2])
        N2 = int(sys.argv[3])
        M2 = int(sys.argv[4])
        if not all(x > 0 for x in [N1, M1, N2, M2]):
            print("Error: All N and M values must be positive integers.")
            return
    except ValueError:
        print("Error: Please provide valid integers for N1, M1, N2, M2.")
        return

    # Phase 1: Expected time for the first collision (5 particles -> 3 particles)
    # Rate lambda_1 = 1. Number of particles n = 5.
    # Formula for E[tau_1] is (1/((n-2)*lambda)) * sum_{1<=i<j<=n-1} d_i * d_j
    # Here, n=5, lambda=1, so pre-factor is 1/3.
    # The gaps are d1=N1, d2=M1, d3=N2, d4=M2.
    
    sum_of_gap_products = N1*M1 + N1*N2 + N1*M2 + M1*N2 + M1*M2 + N2*M2
    
    E_tau1_num = sum_of_gap_products
    E_tau1_den = 3
    
    # Phase 2: Expected time for the second collision (3 particles -> 1 particle)
    # Rate lambda_2 = 2. Number of particles n = 3.
    # The expected time E[tau_2] is given by (1/2) * E[d'_1 * d'_2],
    # where d'_1, d'_2 are the new gaps.
    # E[d'_1 * d'_2] is the sum of products of non-adjacent initial gaps.
    # Non-adjacent pairs are (d1,d3), (d1,d4), (d2,d4).
    
    expected_new_gap_product = N1*N2 + N1*M2 + M1*M2
    
    E_tau2_num = expected_new_gap_product
    E_tau2_den = 2
    
    # Total expected time E[tau] = E[tau_1] + E[tau_2]
    # To add the fractions, we find a common denominator (6)
    
    total_numerator = (E_tau1_num * 2) + (E_tau2_num * 3)
    total_denominator = 6
    
    final_expectation = total_numerator / total_denominator
    
    # Outputting the equation with the final result
    # E[τ] = (1/3) * (N1*M1 + N1*N2 + N1*M2 + M1*N2 + M1*M2 + N2*M2) + 
    #        (1/2) * (N1*N2 + N1*M2 + M1*M2)
    
    term1 = f"(1/3) * ({N1}*{M1} + {N1}*{N2} + {N1}*{M2} + {M1}*{N2} + {M1}*{M2} + {N2}*{M2})"
    term2 = f"(1/2) * ({N1}*{N2} + {N1}*{M2} + {M1}*{M2})"
    
    print(f"The total expected time E[τ] is calculated by the formula:")
    print(f"E[τ] = {term1} + {term2}")
    
    val_term1 = (1/3) * sum_of_gap_products
    val_term2 = (1/2) * expected_new_gap_product
    
    print(f"E[τ] = {val_term1:.4f} + {val_term2:.4f}")
    print(f"E[τ] = {final_expectation:.4f}")

if __name__ == "__main__":
    solve()