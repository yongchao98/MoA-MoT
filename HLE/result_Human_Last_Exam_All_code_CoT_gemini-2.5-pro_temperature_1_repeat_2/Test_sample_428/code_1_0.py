def solve_expectation(N1, M1, N2, M2):
    """
    Calculates the expected time until only one particle remains.

    Args:
        N1, M1, N2, M2: Positive integers defining the initial particle separations.
    """

    # Stage 1: Time from 5 particles to 3. Rate lambda_1 = 1.
    # E[tau_1] = (1/3) * (N1*M1 + N1*N2 + N1*M2 + M1*N2 + M1*M2 + N2*M2)
    e_tau1_num = N1*M1 + N1*N2 + N1*M2 + M1*N2 + M1*M2 + N2*M2
    e_tau1_den = 3
    
    # Stage 2: Time from 3 particles to 1. Rate lambda_2 = 2.
    # This is averaged over all 10 possible 3-particle subsystems.
    # The sum of products of gaps is:
    # 3*N1*M1 + 4*N1*N2 + 3*N1*M2 + 4*M1*N2 + 4*M1*M2 + 3*N2*M2
    e_tau2_num_sum_term = 3*N1*M1 + 4*N1*N2 + 3*N1*M2 + 4*M1*N2 + 4*M1*M2 + 3*N2*M2
    # Denominator is C(5,3) * lambda_2 = 10 * 2 = 20
    e_tau2_den = 20
    
    # Total expectation E[tau] = E[tau_1] + E[tau_2]
    # We combine the expressions by finding a common denominator (60).
    # E[tau] = (20 * e_tau1_num + 3 * e_tau2_num_sum_term) / 60
    
    # Coefficients for each term N_i*M_j after combining
    # N1*M1: 1/3 + 3/20 = (20+9)/60 = 29/60
    # N1*N2: 1/3 + 4/20 = (20+12)/60 = 32/60
    # N1*M2: 1/3 + 3/20 = (20+9)/60 = 29/60
    # M1*N2: 1/3 + 4/20 = (20+12)/60 = 32/60
    # M1*M2: 1/3 + 4/20 = (20+12)/60 = 32/60
    # N2*M2: 1/3 + 3/20 = (20+9)/60 = 29/60

    term1 = N1 * M1
    term2 = N1 * N2
    term3 = N1 * M2
    term4 = M1 * N2
    term5 = M1 * M2
    term6 = N2 * M2
    
    total_expectation = (29/60) * term1 + (32/60) * term2 + \
                        (29/60) * term3 + (32/60) * term4 + \
                        (32/60) * term5 + (29/60) * term6
                        
    print("The expectation E[tau] is given by the formula:")
    print("E[tau] = (29/60)*N1*M1 + (32/60)*N1*N2 + (29/60)*N1*M2 + (32/60)*M1*N2 + (32/60)*M1*M2 + (29/60)*N2*M2")
    print("\nFor the given values:")
    print(f"N1 = {N1}, M1 = {M1}, N2 = {N2}, M2 = {M2}")
    print("\nThe contributions to the formula are:")
    print(f"(29/60) * ({N1} * {M1}) = {(29/60) * term1:.2f}")
    print(f"(32/60) * ({N1} * {N2}) = {(32/60) * term2:.2f}")
    print(f"(29/60) * ({N1} * {M2}) = {(29/60) * term3:.2f}")
    print(f"(32/60) * ({M1} * {N2}) = {(32/60) * term4:.2f}")
    print(f"(32/60) * ({M1} * {M2}) = {(32/60) * term5:.2f}")
    print(f"(29/60) * ({N2} * {M2}) = {(29/60) * term6:.2f}")
    print("\nThe final expectation is:")
    print(f"E[tau] = {total_expectation}")


# Example values for N1, M1, N2, M2
N1_val = 10
M1_val = 5
N2_val = 10
M2_val = 5

solve_expectation(N1_val, M1_val, N2_val, M2_val)
