def solve_random_walk():
    """
    Calculates the expected time until only one particle remains.

    The user should modify the values of N1, M1, N2, M2 as needed.
    """
    # Positive integers defining the initial particle locations.
    # You can change these values to solve for your specific problem.
    N1 = 10
    M1 = 10
    N2 = 10
    M2 = 10

    # Rates for the two phases of the process
    lambda1 = 1.0
    lambda2 = 2.0

    # Sum of products of adjacent initial gaps
    S_adj = N1 * M1 + M1 * N2 + N2 * M2
    
    # Sum of products of non-adjacent initial gaps
    S_non_adj = N1 * N2 + N1 * M2 + M1 * M2
    
    # The formula for the expected total time
    expected_tau = S_non_adj / lambda1 + S_adj / lambda2

    print("Given the initial separations:")
    print(f"N1 = {N1}")
    print(f"M1 = {M1}")
    print(f"N2 = {N2}")
    print(f"M2 = {M2}")
    print("\nThe calculation is based on the formula: E[tau] = (N1*N2 + N1*M2 + M1*M2)/lambda1 + (N1*M1 + M1*N2 + N2*M2)/lambda2")
    
    # Breaking down the calculation for clarity
    term1_val = N1*N2 + N1*M2 + M1*M2
    term2_val = N1*M1 + M1*N2 + N2*M2
    
    print(f"\nE[tau] = ({N1}*{N2} + {N1}*{M2} + {M1}*{M2}) / {lambda1} + ({N1}*{M1} + {M1}*{N2} + {N2}*{M2}) / {lambda2}")
    print(f"E[tau] = ({term1_val}) / {lambda1} + ({term2_val}) / {lambda2}")
    print(f"E[tau] = {term1_val / lambda1} + {term2_val / lambda2}")
    print(f"E[tau] = {expected_tau}")
    print("\nFinal Answer:")
    print(expected_tau)

solve_random_walk()