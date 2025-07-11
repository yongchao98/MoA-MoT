def calculate_expected_time(N1, M1, N2, M2):
    """
    Calculates the expected time until only one particle remains.

    Args:
        N1, M1, N2, M2: Positive integers defining the initial particle positions.
    """
    if not all(isinstance(i, int) and i > 0 for i in [N1, M1, N2, M2]):
        print("Error: N1, M1, N2, and M2 must be positive integers.")
        return

    # Rate for the first phase (5 particles)
    lambda1 = 1
    # Rate for the second phase (3 particles)
    lambda2 = 2

    # Calculate expectation of the first phase duration
    E_tau1_num = N1 * M1 + N2 * M2
    E_tau1_den = 2 * lambda1
    E_tau1 = E_tau1_num / E_tau1_den

    # Calculate expectation of the second phase duration
    E_tau2_num = (N1 + M1) * (N2 + M2)
    E_tau2_den = lambda2
    E_tau2 = E_tau2_num / E_tau2_den

    # Total expected time
    E_tau = E_tau1 + E_tau2

    # Print the step-by-step calculation
    print("The expectation of the time until the second collision, E[tau], can be calculated in two parts.")
    print("Part 1: Expected time for the first collision (5 -> 3 particles), E[tau1]")
    print(f"E[tau1] = (N1 * M1 + N2 * M2) / (2 * lambda1)")
    print(f"E[tau1] = ({N1} * {M1} + {N2} * {M2}) / (2 * {lambda1}) = {E_tau1}")
    print("\nPart 2: Expected time for the second collision (3 -> 1 particle), E[tau2]")
    print(f"E[tau2] = ((N1 + M1) * (N2 + M2)) / lambda2")
    print(f"E[tau2] = (({N1} + {M1}) * ({N2} + {M2})) / {lambda2} = {E_tau2}")
    print("\nTotal Expected Time, E[tau] = E[tau1] + E[tau2]")
    print(f"E[tau] = {E_tau1} + {E_tau2} = {E_tau}")

# Example values for N1, M1, N2, M2
N1 = 10
M1 = 5
N2 = 8
M2 = 4

calculate_expected_time(N1, M1, N2, M2)
<<<333.0>>>