import numpy as np

def solve_whitening_filter():
    """
    Calculates and prints the expression for the whitening filter W(D)
    based on the corrected form of the problem.
    """
    
    # Based on the analysis that the original problem for q_k is ill-posed,
    # we assume the intended channel's causal minimum-phase spectral factor is:
    # G(D) = (1 - 1/3 D) / (1 - 1/2 D).
    # The required whitening filter W(D) that makes the resulting channel causal
    # is W(D) = 1 / G(D^{-1}).
    
    # Numerator coefficients for W(D) in powers of D^{-1}
    num_coeff_1 = 1.0
    num_coeff_2 = -1.0 / 2.0
    
    # Denominator coefficients for W(D) in powers of D^{-1}
    den_coeff_1 = 1.0
    den_coeff_2 = -1.0 / 3.0
    
    print("The appropriate whitening filter W(D) that ensures the resulting communication channel is causal is given by the equation:")
    print("")
    
    # Printing the full equation with signs
    print(f"      {num_coeff_1} + ({num_coeff_2}) D^(-1)")
    print(f"W(D) = -----------------------------")
    print(f"      {den_coeff_1} + ({den_coeff_2}) D^(-1)")
    
    print("\nThe individual numbers in the final equation are:")
    print(f"Numerator Constant Term: {num_coeff_1}")
    print(f"Numerator D^-1 Coefficient: {num_coeff_2}")
    print(f"Denominator Constant Term: {den_coeff_1}")
    print(f"Denominator D^-1 Coefficient: {den_coeff_2}")

solve_whitening_filter()
