import sys

def solve_particle_annihilation():
    """
    Calculates the expected time until only one particle remains.
    
    The user should provide four positive integers N1, M1, N2, M2 as command-line arguments.
    If not provided, example values will be used.
    """
    try:
        # Use command-line arguments if provided
        if len(sys.argv) == 5:
            N1 = int(sys.argv[1])
            M1 = int(sys.argv[2])
            N2 = int(sys.argv[3])
            M2 = int(sys.argv[4])
            print(f"Using provided values: N1={N1}, M1={M1}, N2={N2}, M2={M2}")
        else:
            # Example values
            N1, M1, N2, M2 = 2, 3, 4, 5
            print(f"Using example values: N1={N1}, M1={M1}, N2={N2}, M2={M2}")
            print("You can provide your own values as command-line arguments.")
            print("Usage: python your_script_name.py N1 M1 N2 M2\n")

        if not all(isinstance(i, int) and i > 0 for i in [N1, M1, N2, M2]):
            raise ValueError("All inputs must be positive integers.")

    except (ValueError, IndexError):
        print("Error: Please provide four positive integers for N1, M1, N2, M2.")
        print("Usage: python your_script_name.py N1 M1 N2 M2")
        return

    # Phase 1: 5 particles -> 3 particles, rate lambda=1
    # E[tau_1] = (1/2) * (N1*M1 + M1*N2 + N2*M2)
    expected_tau1 = 0.5 * (N1 * M1 + M1 * N2 + N2 * M2)

    # Phase 2: 3 particles -> 1 particle, rate lambda=2
    # E[tau_2] = (1/4) * (N1+M1)*(N2+M2)
    expected_tau2 = 0.25 * (N1 + M1) * (N2 + M2)

    # Total expected time E[tau] = E[tau_1] + E[tau_2]
    total_expected_time = expected_tau1 + expected_tau2

    # Output the breakdown of the calculation
    print("The calculation is based on the formula E[tau] = E[tau1] + E[tau2]")
    print("E[tau] = (1/2)*(N1*M1 + M1*N2 + N2*M2) + (1/4)*(N1+M1)*(N2+M2)")
    print("\nStep-by-step calculation:")
    
    # Print E[tau1] calculation
    print(f"E[tau1] = (1/2)*({N1}*{M1} + {M1}*{N2} + {N2}*{M2})")
    print(f"E[tau1] = 0.5 * ({N1*M1} + {M1*N2} + {N2*M2})")
    print(f"E[tau1] = 0.5 * {N1*M1 + M1*N2 + N2*M2}")
    print(f"E[tau1] = {expected_tau1}")

    # Print E[tau2] calculation
    print(f"\nE[tau2] = (1/4)*({N1}+{M1})*({N2}+{M2})")
    print(f"E[tau2] = 0.25 * ({N1+M1}) * ({N2+M2})")
    print(f"E[tau2] = 0.25 * { (N1+M1)*(N2+M2) }")
    print(f"E[tau2] = {expected_tau2}")

    # Print the final result
    print(f"\nTotal Expected Time E[tau] = {expected_tau1} + {expected_tau2}")
    print(f"E[tau] = {total_expected_time}")
    
    # Final answer in the required format
    # Note: The 'answer' format is specified for a single value.
    # The full calculation is more informative, but if a single number is required,
    # it would be total_expected_time.
    # print(f"\n<<<{total_expected_time}>>>")

if __name__ == "__main__":
    solve_particle_annihilation()