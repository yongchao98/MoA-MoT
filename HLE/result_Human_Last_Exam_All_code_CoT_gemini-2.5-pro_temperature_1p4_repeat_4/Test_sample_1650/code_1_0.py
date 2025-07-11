import sympy

def display_overlap_integral_2s():
    """
    This function presents the derived analytical expression for the overlap
    integral of two 2s orbitals in the H2+ ion.
    """

    # Define the symbols for clarity in the equation.
    # R represents the internuclear distance.
    # zeta (ζ) represents the effective nuclear charge.
    R = sympy.Symbol('R')
    zeta = sympy.Symbol('ζ')
    
    # For compactness, the formula is often written using rho (ρ) = ζ * R.
    rho = zeta * R
    
    # The derived analytical expression for the 2s-2s overlap integral (S) is:
    # S(ρ) = exp(-ρ/2) * (1 + ρ/2 + ρ²/12 - ρ⁴/240)
    # We will now construct and print this formula in terms of R and ζ.
    
    term0 = 1
    term1_num = zeta * R
    term1_den = 2
    term2_num = (zeta * R)**2
    term2_den = 12
    term3_num = (zeta * R)**4
    term3_den = 240
    
    print("The overlap integral for two 2s orbitals in H₂⁺ is given by S(R, ζ).")
    print("The analytical expression is derived by evaluating the integral ∫ψ₂s(A)*ψ₂s(B)dτ in elliptical coordinates.")
    print("\nThe final formula is:")
    
    # Output the final equation, showing each number as requested.
    equation_str = (
        f"S(R, ζ) = exp(-({zeta} * {R}) / {term1_den}) * "
        f"({term0} + ({term1_num}) / {term1_den} + ({term2_num}) / {term2_den} - ({term3_num}) / {term3_den})"
    )
    
    print(equation_str)

if __name__ == '__main__':
    display_overlap_integral_2s()