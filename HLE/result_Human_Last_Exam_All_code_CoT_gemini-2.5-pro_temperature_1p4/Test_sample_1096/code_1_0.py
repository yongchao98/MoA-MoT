import math

def solve_and_print_formulas():
    """
    This function calculates the coefficients for the optimal step sizes
    gamma_1 and gamma_2 and prints the final formulas.
    """
    
    # The coefficients are derived from the roots of the Chebyshev polynomial
    sqrt2 = math.sqrt(2)
    
    # Coefficients for the numerator of gamma_1: (4-2*sqrt(2))*kappa + (4+2*sqrt(2))
    a1 = 4 - 2 * sqrt2
    b1 = 4 + 2 * sqrt2
    
    # Coefficients for the numerator of gamma_2: (4+2*sqrt(2))*kappa + (4-2*sqrt(2))
    a2 = 4 + 2 * sqrt2
    b2 = 4 - 2 * sqrt2
    
    # The denominator is kappa^2 + 6*kappa + 1 for both
    c = 6
    d = 1
    
    print("The optimal choice for the pair (gamma_1, gamma_2) is given by the following formulas, where kappa is the condition number M/mu:")
    print("-" * 80)
    
    # Print the formula for gamma_1
    print(f"gamma_1 = ({a1:.8f} * kappa + {b1:.8f}) / (kappa^2 + {c} * kappa + {d})")
    
    # Print the formula for gamma_2
    print(f"gamma_2 = ({a2:.8f} * kappa + {b2:.8f}) / (kappa^2 + {c} * kappa + {d})")
    print("-" * 80)
    print("Note: The order of gamma_1 and gamma_2 can be swapped.")

# Execute the function to display the answer
solve_and_print_formulas()