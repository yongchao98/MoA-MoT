import sys

def solve_particle_problem(N1, M1, N2, M2):
    """
    Calculates the expected time for the second collision of annihilating random walks.

    Args:
        N1 (int): First gap size.
        M1 (int): Second gap size.
        N2 (int): Third gap size.
        M2 (int): Fourth gap size.
    """
    # Check if inputs are positive integers
    if not all(isinstance(i, int) and i > 0 for i in [N1, M1, N2, M2]):
        print("Error: N1, M1, N2, M2 must be positive integers.", file=sys.stderr)
        return

    # Phase 1: 5 particles -> 3 particles, rate lambda_1 = 1
    # The expected time is (1/3) * sum of all pairwise products of initial gaps.
    sum_of_gap_products = N1*M1 + N1*N2 + N1*M2 + M1*N2 + M1*M2 + N2*M2
    e_tau1 = (1/3) * sum_of_gap_products

    # Phase 2: 3 particles -> 1 particle, rate lambda_2 = 2
    # The expected time for the second phase is (1/lambda_2) * E[product of new gaps].
    # E[product of new gaps] is the sum of products of non-adjacent initial gaps.
    lambda2 = 2
    sum_of_non_adjacent_gap_products = N1*N2 + N1*M2 + M1*M2
    e_tau2 = sum_of_non_adjacent_gap_products / lambda2

    # Total expected time
    total_expectation = e_tau1 + e_tau2
    
    # We express the final answer as a sum of two parts, E[tau1] and E[tau2].
    # E[tau] = 1/3 * (N1*M1 + N1*N2 + N1*M2 + M1*N2 + M1*M2 + N2*M2) + 1/2 * (N1*N2 + N1*M2 + M1*M2)
    # The output will show each term in the final equation.
    
    print(f"Let N1 = {N1}, M1 = {M1}, N2 = {N2}, M2 = {M2}.")
    print("The expectation of tau is E[tau] = E[tau1] + E[tau2].")
    print("\nE[tau1] = 1/3 * (N1*M1 + N1*N2 + N1*M2 + M1*N2 + M1*M2 + N2*M2)")
    print(f"E[tau1] = 1/3 * ({N1}*{M1} + {N1}*{N2} + {N1}*{M2} + {M1}*{N2} + {M1}*{M2} + {N2}*{M2})")
    print(f"E[tau1] = 1/3 * ({N1*M1} + {N1*N2} + {N1*M2} + {M1*N2} + {M1*M2} + {N2*M2})")
    print(f"E[tau1] = 1/3 * {sum_of_gap_products} = {e_tau1:.4f}")
    
    print("\nE[tau2] = 1/2 * (N1*N2 + N1*M2 + M1*M2)")
    print(f"E[tau2] = 1/2 * ({N1}*{N2} + {N1}*{M2} + {M1}*{M2})")
    print(f"E[tau2] = 1/2 * ({N1*N2} + {N1*M2} + {M1*M2})")
    print(f"E[tau2] = 1/2 * {sum_of_non_adjacent_gap_products} = {e_tau2:.4f}")
    
    print(f"\nTotal Expectation E[tau] = {e_tau1:.4f} + {e_tau2:.4f} = {total_expectation:.4f}")


# Example usage with some fixed positive integers, e.g., N1=1, M1=2, N2=3, M2=4
# You can change these values to test with other numbers.
if __name__ == '__main__':
    # You can provide N1, M1, N2, M2 as command-line arguments.
    # e.g., python your_script.py 1 2 3 4
    if len(sys.argv) == 5:
        try:
            N1, M1, N2, M2 = map(int, sys.argv[1:])
            solve_particle_problem(N1, M1, N2, M2)
        except ValueError:
            print("Please provide four integer arguments for N1, M1, N2, M2.", file=sys.stderr)
    else:
        # Default example if no arguments are given.
        print("Using default example values: N1=1, M1=2, N2=3, M2=4\n")
        solve_particle_problem(1, 2, 3, 4)
