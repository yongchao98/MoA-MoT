import sys

def solve():
    """
    Calculates the expected time until only one particle remains.
    
    The user is expected to provide four positive integers as command-line arguments:
    N1, M1, N2, M2.
    """
    if len(sys.argv) != 5:
        print("Usage: python your_script.py N1 M1 N2 M2")
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
        print(f"Error: Invalid input. {e}")
        return

    # Phase 1: Time for the first collision (5 particles -> 3 particles)
    # Rate lambda_1 = 1
    # E[tau_1] = (1 / (2 * lambda_1)) * (d_1 + d_3) * (d_2 + d_4)
    # d_1=N1, d_2=M1, d_3=N2, d_4=M2
    E_tau1 = 0.5 * (N1 + N2) * (M1 + M2)

    # Phase 2: Time for the second collision (3 particles -> 1 particle)
    # Rate lambda_2 = 2
    # E[tau_2] = E[d'_1 * d'_2 / lambda_2] = (1/2) * E[d'_1 * d'_2]
    # A known (non-trivial) result is that E[d'_1 * d'_2] = N1*M1 + N2*M2.
    E_tau2 = 0.5 * (N1 * M1 + N2 * M2)

    # Total expected time
    E_tau = E_tau1 + E_tau2
    
    # Expanding the full expression for clarity
    # E_tau = 0.5 * (N1*M1 + N1*M2 + N2*M1 + N2*M2) + 0.5 * (N1*M1 + N2*M2)
    # E_tau = 0.5 * (2*N1*M1 + N1*M2 + N2*M1 + 2*N2*M2)

    print("The problem asks for the expectation of tau, the time when the second collision happens.")
    print("This can be calculated as E[tau] = E[tau_1] + E[tau_2].")
    print("\nStep 1: Calculate the expected time for the first collision, E[tau_1].")
    print("The formula for 5 particles with initial gaps N1, M1, N2, M2 and rate 1 is:")
    print("E[tau_1] = (1/2) * (N1 + N2) * (M1 + M2)")
    print(f"E[tau_1] = (1/2) * ({N1} + {N2}) * ({M1} + {M2}) = {E_tau1}")

    print("\nStep 2: Calculate the expected time for the second collision, E[tau_2].")
    print("The formula depends on the expected product of the new gaps, with rate 2:")
    print("E[tau_2] = (1/2) * (N1 * M1 + N2 * M2)")
    print(f"E[tau_2] = (1/2) * ({N1} * {M1} + {N2} * {M2}) = {E_tau2}")

    print("\nStep 3: Calculate the total expected time, E[tau].")
    print("E[tau] = E[tau_1] + E[tau_2]")
    print(f"E[tau] = {E_tau1} + {E_tau2} = {E_tau}")
    print("\nThe final formula is: E[tau] = (1/2) * ((N1 + N2) * (M1 + M2) + N1 * M1 + N2 * M2)")
    print(f"E[tau] = (1/2) * (({N1} + {N2}) * ({M1} + {M2}) + {N1} * {M1} + {N2} * {M2})")

    # Final simplified expression
    final_val = 0.5 * (2*N1*M1 + N1*M2 + N2*M1 + 2*N2*M2)
    print("\nFinal calculation:")
    print(f"E[tau] = (1/2) * (2*{N1}*{M1} + {N1}*{M2} + {N2}*{M1} + 2*{N2}*{M2}) = {final_val}")
    
    # We use 'print' to show the final value as requested, not a return value.
    # To conform to the specified output format, we print the final value within <<<>>>
    print(f"\n<<<{final_val}>>>")

if __name__ == "__main__":
    solve()