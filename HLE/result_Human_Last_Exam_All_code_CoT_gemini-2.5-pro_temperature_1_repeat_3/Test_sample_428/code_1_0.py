import sys

def solve_annihilation_time(N1, M1, N2, M2):
    """
    Calculates the expected time until only one particle remains.
    
    The formula is derived from the theory of annihilating random walks.
    Let E1 be the expected time for the first collision (5 -> 3 particles, rate 1).
    Let E2 be the expected time for the second collision (3 -> 1 particle, rate 2).
    E1 = N1*M1 + M1*N2 + N2*M2
    E2 = (N1+M1)*(N2+M2) / 2
    Total expected time E = E1 + E2.
    
    Args:
        N1, M1, N2, M2: Positive integers defining the initial particle positions.
    """
    if not all(isinstance(i, int) and i > 0 for i in [N1, M1, N2, M2]):
        print("Error: N1, M1, N2, and M2 must be positive integers.")
        return

    # Phase 1: 5 particles -> 3 particles (rate lambda=1)
    e1 = N1 * M1 + M1 * N2 + N2 * M2
    
    # Phase 2: 3 particles -> 1 particle (rate lambda=2)
    e2 = (N1 + M1) * (N2 + M2) / 2
    
    # Total expected time
    total_expectation = e1 + e2
    
    print("The problem is to find the expectation of tau, where tau is the time of the second collision.")
    print("The expectation of tau is given by the formula E[tau] = E[tau_1] + E[tau - tau_1]")
    print("\nPhase 1: Time for the first collision (5 -> 3 particles, rate=1)")
    print(f"E[tau_1] = N1*M1 + M1*N2 + N2*M2")
    print(f"E[tau_1] = {N1}*{M1} + {M1}*{N2} + {N2}*{M2} = {e1}")

    print("\nPhase 2: Time for the second collision (3 -> 1 particle, rate=2)")
    print(f"E[tau - tau_1] = (N1+M1)*(N2+M2) / 2")
    print(f"E[tau - tau_1] = ({N1}+{M1})*({N2}+{M2}) / 2 = {e2}")

    print("\nTotal expected time:")
    print(f"E[tau] = ({N1}*{M1} + {M1}*{N2} + {N2}*{M2}) + (({N1}+{M1})*({N2}+{M2}))/2")
    print(f"E[tau] = {e1} + {e2} = {total_expectation}")

    return total_expectation

if __name__ == '__main__':
    # Example usage with some placeholder values.
    # The user should replace these with their specific values for N1, M1, N2, M2.
    try:
        # Check if command-line arguments are provided
        if len(sys.argv) == 5:
            N1 = int(sys.argv[1])
            M1 = int(sys.argv[2])
            N2 = int(sys.argv[3])
            M2 = int(sys.argv[4])
        else:
            # Default values if no arguments are given
            print("Usage: python your_script.py N1 M1 N2 M2")
            print("Using default example values N1=1, M1=1, N2=1, M2=1\n")
            N1, M1, N2, M2 = 1, 1, 1, 1
        
        result = solve_annihilation_time(N1, M1, N2, M2)
    except (ValueError, IndexError):
        print("Invalid input. Please provide four positive integers for N1, M1, N2, and M2.")
