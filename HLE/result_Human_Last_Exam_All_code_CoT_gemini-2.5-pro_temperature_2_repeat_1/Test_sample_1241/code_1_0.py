import numpy as np

def solve_kolmogorov_steady_state():
    """
    Solves the Kolmogorov-Chapman system for steady-state probabilities and calculates P0 + P1.
    """
    # Given transition rates
    lambda01 = 0.019
    lambda10 = 0.65
    lambda12 = 0.4
    lambda21 = 0.392
    lambda23 = 0.008
    lambda31 = 0.008

    # At steady state (t -> +inf), derivatives P_i'(t) are 0.
    # The system becomes a set of linear equations for probabilities p_i = P_i(+inf).
    #
    # From P_0'(t) = 0:  lambda01 * p0 = lambda10 * p1  => p0 = (lambda10 / lambda01) * p1
    # From P_2'(t) = 0:  (lambda21 + lambda23) * p2 = lambda12 * p1 => p2 = (lambda12 / (lambda21 + lambda23)) * p1
    # From P_3'(t) = 0:  lambda31 * p3 = lambda23 * p2 => p3 = (lambda23 / lambda31) * p2
    #
    # We express p0, p2, and p3 as multiples of p1. Let's find the coefficients.
    
    # Coefficient for p0 in terms of p1
    coeff_p0 = lambda10 / lambda01
    
    # Coefficient for p2 in terms of p1
    coeff_p2 = lambda12 / (lambda21 + lambda23)
    
    # Coefficient for p3 in terms of p1 (by first relating p3 to p2, then p2 to p1)
    coeff_p3_of_p2 = lambda23 / lambda31
    coeff_p3 = coeff_p3_of_p2 * coeff_p2
    
    # Now use the normalization condition: p0 + p1 + p2 + p3 = 1
    # (coeff_p0 * p1) + p1 + (coeff_p2 * p1) + (coeff_p3 * p1) = 1
    # p1 * (coeff_p0 + 1 + coeff_p2 + coeff_p3) = 1
    
    # Solve for p1
    total_coeffs = coeff_p0 + 1 + coeff_p2 + coeff_p3
    p1 = 1 / total_coeffs
    
    # Calculate p0
    p0 = coeff_p0 * p1
    
    # Calculate the required sum P0(inf) + P1(inf)
    sum_p0_p1 = p0 + p1
    
    # Print the result in the required format
    print(f"P_0(+\infty) + P_1(+\infty) = {p0} + {p1} = {sum_p0_p1}")

# Execute the function to print the result
solve_kolmogorov_steady_state()