import numpy as np

def solve_markov_chain():
    """
    Solves for the stationary probabilities of the given Kolmogorov-Chepmen system
    and calculates the sum P0 + P1.
    """
    # Given transition rates
    lambda_01 = 0.019
    lambda_10 = 0.65
    lambda_12 = 0.4
    lambda_21 = 0.392
    lambda_23 = 0.008
    lambda_31 = 0.008

    # At steady state (t -> infinity), P_i'(t) = 0.
    # The system of differential equations becomes a system of linear equations.
    # We will express P0, P2, and P3 in terms of P1.

    # From -l_01*P0 + l_10*P1 = 0
    # P0 = (l_10 / l_01) * P1
    coeff_P0_wrt_P1 = lambda_10 / lambda_01

    # From l_12*P1 - (l_21 + l_23)*P2 = 0
    # P2 = (l_12 / (l_21 + l_23)) * P1
    coeff_P2_wrt_P1 = lambda_12 / (lambda_21 + lambda_23)

    # From l_23*P2 - l_31*P3 = 0
    # P3 = (l_23 / l_31) * P2
    # Substitute P2 in terms of P1
    # P3 = (l_23 / l_31) * coeff_P2_wrt_P1 * P1
    coeff_P3_wrt_P1 = (lambda_23 / lambda_31) * coeff_P2_wrt_P1

    # The normalization condition is P0 + P1 + P2 + P3 = 1.
    # Substitute the expressions in terms of P1:
    # (coeff_P0_wrt_P1 * P1) + P1 + (coeff_P2_wrt_P1 * P1) + (coeff_P3_wrt_P1 * P1) = 1
    # P1 * (coeff_P0_wrt_P1 + 1 + coeff_P2_wrt_P1 + coeff_P3_wrt_P1) = 1
    
    total_coeff_P1 = coeff_P0_wrt_P1 + 1 + coeff_P2_wrt_P1 + coeff_P3_wrt_P1
    
    # Solve for P1
    P1 = 1 / total_coeff_P1
    
    # Calculate P0
    P0 = coeff_P0_wrt_P1 * P1

    # The desired result is the sum of P0 and P1
    result = P0 + P1

    # Print the final equation with the calculated values
    print(f"The stationary probabilities are:")
    print(f"P0(inf) = {P0}")
    print(f"P1(inf) = {P1}")
    print("\nThe final equation is:")
    print(f"P0(inf) + P1(inf) = {P0} + {P1} = {result}")

solve_markov_chain()