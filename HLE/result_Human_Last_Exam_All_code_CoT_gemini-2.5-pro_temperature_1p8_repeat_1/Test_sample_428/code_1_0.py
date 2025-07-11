import sys

def solve():
    """
    Calculates the expected time until only one particle remains.
    
    The user is expected to provide four positive integer values for
    N1, M1, N2, M2 as command-line arguments.
    """
    if len(sys.argv) != 5:
        print("Usage: python your_script_name.py N1 M1 N2 M2")
        print("Please provide four positive integers for N1, M1, N2, and M2.")
        return

    try:
        N1 = int(sys.argv[1])
        M1 = int(sys.argv[2])
        N2 = int(sys.argv[3])
        M2 = int(sys.argv[4])
        
        if not all(x > 0 for x in [N1, M1, N2, M2]):
            raise ValueError("Inputs must be positive integers.")

    except ValueError as e:
        print(f"Error: {e}")
        print("Please ensure you provide four positive integers as command-line arguments.")
        return

    # Expected time for the first collision (5 particles -> 3 particles)
    # Rate lambda_1 = 1
    # E[tau_1] is approximated by the sum of expected collision times of adjacent 3-particle subsystems.
    E_tau1_term1 = N1 * M1
    E_tau1_term2 = M1 * N2
    E_tau1_term3 = N2 * M2
    E_tau1 = E_tau1_term1 + E_tau1_term2 + E_tau1_term3
    
    # Expected time for the second collision (3 particles -> 1 particle)
    # Rate lambda_2 = 2
    # This is modeled as the collision time for a 3-particle system composed of
    # the outer particles (1 and 5) and the central particle (3).
    # The gaps are (N1+M1) and (N2+M2).
    E_tau2 = (N1 + M1) * (N2 + M2) / 2.0

    # Total expected time
    E_tau = E_tau1 + E_tau2
    
    # Outputting the final equation step-by-step
    print(f"Given N1 = {N1}, M1 = {M1}, N2 = {N2}, M2 = {M2}:")
    print("\nStep 1: Calculate the expected time for the first collision (E[tau_1])")
    print(f"E[tau_1] = N1*M1 + M1*N2 + N2*M2")
    print(f"E[tau_1] = {N1}*{M1} + {M1}*{N2} + {N2}*{M2}")
    print(f"E[tau_1] = {E_tau1_term1} + {E_tau1_term2} + {E_tau1_term3} = {E_tau1}")
    
    print("\nStep 2: Calculate the expected time for the second collision (E[tau_2])")
    print(f"E[tau_2] = (N1 + M1) * (N2 + M2) / 2")
    print(f"E[tau_2] = ({N1} + {M1}) * ({N2} + {M2}) / 2")
    print(f"E[tau_2] = {N1+M1} * {N2+M2} / 2 = {E_tau2}")

    print("\nStep 3: Calculate the total expected time E[tau]")
    print(f"E[tau] = E[tau_1] + E[tau_2]")
    print(f"E[tau] = {E_tau1} + {E_tau2}")
    print(f"The final expectation of tau is: {E_tau}")


if __name__ == "__main__":
    solve()
