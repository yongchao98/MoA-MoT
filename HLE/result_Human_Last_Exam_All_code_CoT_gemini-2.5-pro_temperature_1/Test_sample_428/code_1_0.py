import sys

def solve():
    """
    Calculates the expected time for the second collision of annihilating random walks.
    """
    try:
        # Read inputs for N1, M1, N2, M2 from the command line
        if len(sys.argv) != 5:
            print("Usage: python your_script.py N1 M1 N2 M2")
            print("Example: python your_script.py 10 20 15 25")
            return

        N1 = int(sys.argv[1])
        M1 = int(sys.argv[2])
        N2 = int(sys.argv[3])
        M2 = int(sys.argv[4])

        if N1 <= 0 or M1 <= 0 or N2 <= 0 or M2 <= 0:
            print("Error: N1, M1, N2, M2 must be positive integers.")
            return

    except (ValueError, IndexError):
        # Default values if command line arguments are not provided or invalid
        print("Using default values N1=2, M1=2, N2=3, M2=3 as an example.")
        N1, M1, N2, M2 = 2, 2, 3, 3

    # The rate for the first phase (5 to 3 particles) is lambda_1 = 1
    # The rate for the second phase (3 to 1 particle) is lambda_2 = 2
    
    # Expected time for the first collision (5 -> 3 particles)
    # E[tau_1] = (N1*M1 + N2*M2) / (2 * lambda_1)
    e_tau1_num = N1 * M1 + N2 * M2
    e_tau1_den = 2
    e_tau1 = e_tau1_num / e_tau1_den
    
    # Expected time for the second collision (3 -> 1 particle)
    # E[tau - tau_1] = E[G1*G2] / lambda_2
    # E[G1*G2] = 0.5 * ((N1+M1)*(N2+M2) + M1*N2)
    # E[tau - tau_1] = (0.5 * ((N1+M1)*(N2+M2) + M1*N2)) / 2
    # E[tau - tau_1] = ((N1+M1)*(N2+M2) + M1*N2) / 4
    e_tau2_num = (N1 + M1) * (N2 + M2) + M1 * N2
    e_tau2_den = 4
    e_tau2 = e_tau2_num / e_tau2_den
    
    # Total expected time
    total_expected_time = e_tau1 + e_tau2

    # Print the formula with the specific numbers
    print("The formula for the total expected time E[tau] is:")
    print("E[tau] = (N1*M1 + N2*M2)/2 + ((N1+M1)*(N2+M2) + M1*N2)/4")
    print("\nPlugging in the given values:")
    
    # Print the equation with substituted values
    part1_str = f"({N1}*{M1} + {N2}*{M2})/{e_tau1_den}"
    part2_str = f"(({N1}+{M1})*({N2}+{M2}) + {M1}*{N2})/{e_tau2_den}"
    
    print(f"E[tau] = {part1_str} + {part2_str}")
    print(f"E[tau] = {e_tau1_num}/{e_tau1_den} + {e_tau2_num}/{e_tau2_den}")
    print(f"E[tau] = {e_tau1} + {e_tau2}")
    
    # Print the final result
    print(f"\nThe final computed value is:")
    print(f"{total_expected_time}")

if __name__ == "__main__":
    solve()