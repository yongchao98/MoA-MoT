import sympy

def solve_whitening_filter():
    """
    This function determines and prints the expression for the whitening filter W(D).
    
    The solution is based on the assumption that the intended causal spectral factor
    of the channel autocorrelation Q(D) is H_c(D) = (1+D)/(1-D/3), due to the
    coefficients q_1 and q_2 matching those derived from this factor.
    
    The whitening filter W(D) is then calculated as 1/H_c(D^{-1}).
    """
    
    # Define the D-transform variable
    D = sympy.Symbol('D')
    
    # Assumed causal spectral factor H_c(D)
    # H_c = (1 + D) / (1 - D/3)
    
    # Corresponding H_c(D^{-1})
    H_c_invD = (1 + D**-1) / (1 - (D**-1)/3)
    
    # Whitening filter W(D) = 1 / H_c(D^{-1})
    W_D = 1 / H_c_invD
    
    # Simplify the expression for W(D)
    W_D_simplified = sympy.simplify(W_D)
    
    # Extract numerator and denominator for printing
    num, den = W_D_simplified.as_numer_denom()
    
    # Get coefficients for the final equation
    num_const = sympy.Poly(num, D**-1).coeffs()[1]
    num_d_inv_coeff = sympy.Poly(num, D**-1).coeffs()[0]
    
    den_const = sympy.Poly(den, D**-1).coeffs()[1]
    den_d_inv_coeff = sympy.Poly(den, D**-1).coeffs()[0]

    print("Based on analysis, the appropriate whitening filter W(D) is:")
    
    # Format the output string for the equation
    # Using sympy.pretty_print for a cleaner mathematical look
    print("\n       1 - (1/3)*D⁻¹")
    print("W(D) = —————————————")
    print("         1 + D⁻¹\n")

    print("The numbers in the final equation are:")
    print(f"Numerator constant term: {num_const}")
    print(f"Numerator D⁻¹ term coefficient: {num_d_inv_coeff}")
    print(f"Denominator constant term: {den_const}")
    print(f"Denominator D⁻¹ term coefficient: {den_d_inv_coeff}")

solve_whitening_filter()
<<<W(D) = (1 - 1/3 * D**-1) / (1 + D**-1)>>>