import numpy as np

def solve_stationary_probabilities():
    """
    Calculates the stationary probabilities P_i(inf) and their sum P0+P1
    based on the corrected Kolmogorov-Chepmen system.
    """
    # Given transition rates
    lambda_01 = 0.019
    lambda_10 = 0.65
    lambda_12 = 0.4
    lambda_21 = 0.392
    lambda_23 = 0.008
    lambda_31 = 0.008

    # In the steady state (t -> +inf), P_i'(t) = 0. Let p_i = P_i(+inf).
    # From the relations derived in the explanation:
    # p1 = p2 = p3
    # p0 = (lambda_10 / lambda_01) * p1
    
    # Using the normalization condition p0 + p1 + p2 + p3 = 1:
    # p1 * (lambda_10 / lambda_01 + 3) = 1
    
    p1 = 1 / (lambda_10 / lambda_01 + 3)
    
    # Calculate p0
    p0 = (lambda_10 / lambda_01) * p1
    
    # The required sum is p0 + p1
    sum_p0_p1 = p0 + p1
    
    print("The final calculation is:")
    # The requirement is to output each number in the final equation.
    print(f"P0(inf) + P1(inf) = {p0} + {p1} = {sum_p0_p1}")

solve_stationary_probabilities()