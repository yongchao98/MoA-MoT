import numpy as np

def solve_steady_state_probabilities():
    """
    Solves the Kolmogorov-Chepmen system for steady-state probabilities
    and calculates P0(inf) + P1(inf).
    """
    # Given lambda values
    lambda_01 = 0.019
    lambda_10 = 0.65
    lambda_12 = 0.4
    lambda_21 = 0.392
    lambda_23 = 0.008
    lambda_31 = 0.008

    # At steady state, P_i'(t) = 0. We get a system of linear equations.
    # We will express p0, p2, and p3 in terms of p1.

    # From 0 = -lambda_01*p0 + lambda_10*p1  =>  p0 = (lambda_10 / lambda_01) * p1
    coeff_p0_in_p1 = lambda_10 / lambda_01

    # From 0 = lambda_12*p1 - (lambda_21 + lambda_23)*p2  =>  p2 = (lambda_12 / (lambda_21 + lambda_23)) * p1
    coeff_p2_in_p1 = lambda_12 / (lambda_21 + lambda_23)

    # From 0 = lambda_23*p2 - lambda_31*p3  =>  p3 = (lambda_23 / lambda_31) * p2
    # Substitute p2 in terms of p1: p3 = (lambda_23 / lambda_31) * coeff_p2_in_p1 * p1
    coeff_p3_in_p1 = (lambda_23 / lambda_31) * coeff_p2_in_p1

    # Now use the normalization equation: p0 + p1 + p2 + p3 = 1
    # (coeff_p0_in_p1 * p1) + p1 + (coeff_p2_in_p1 * p1) + (coeff_p3_in_p1 * p1) = 1
    # p1 * (coeff_p0_in_p1 + 1 + coeff_p2_in_p1 + coeff_p3_in_p1) = 1
    
    total_coeff = coeff_p0_in_p1 + 1 + coeff_p2_in_p1 + coeff_p3_in_p1
    
    # Solve for p1
    p1 = 1 / total_coeff
    
    # Solve for p0
    p0 = coeff_p0_in_p1 * p1
    
    # The required sum is p0 + p1
    result = p0 + p1
    
    print(f"The steady-state probabilities are:")
    print(f"P0(inf) = {p0}")
    print(f"P1(inf) = {p1}")
    print("\nThe final equation is:")
    print(f"{p0} + {p1} = {result}")

solve_steady_state_probabilities()