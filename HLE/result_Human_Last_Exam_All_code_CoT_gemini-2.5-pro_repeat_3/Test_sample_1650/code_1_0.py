import sympy

def print_overlap_integral_2s_2s():
    """
    This function prints the derived analytical expression for the overlap
    integral (S) of two 2s orbitals.

    The expression is given in terms of rho, where rho = zeta * R.
    - S: The overlap integral
    - zeta (ζ): The effective nuclear charge
    - R: The internuclear distance
    """

    # Using sympy for pretty printing the mathematical expression
    zeta, R = sympy.symbols('ζ R')
    rho = sympy.Symbol('ρ') # Represents zeta * R

    # The coefficients of the polynomial part of the expression
    c0 = 1
    c1 = sympy.Rational(1, 2)
    c2 = sympy.Rational(1, 12)
    c4 = sympy.Rational(1, 240)

    # Constructing the expression S(rho)
    # The derived formula is S(ρ) = exp(-ρ/2) * (1 + ρ/2 + ρ²/12 + ρ⁴/240)
    S_rho = sympy.exp(-rho/2) * (c0 + c1*rho + c2*rho**2 + c4*rho**4)

    # Explanation text
    print("The analytical expression for the overlap integral S of two 2s orbitals is:")
    print("S(R) = exp(-ζR/2) * [c₀ + c₁(ζR) + c₂(ζR)² + c₃(ζR)³ + c₄(ζR)⁴]")
    print("\nWhere ρ = ζR:")
    print(f"S(ρ) = exp(-ρ/2) * [{c0} + ({c1})ρ + ({c2})ρ² + ({c4})ρ⁴]")
    
    print("\nIn formatted mathematical notation:")
    # Using sympy's pretty print feature
    sympy.pprint(S_rho, use_unicode=True)
    
    print("\nFinal equation with explicit numbers:")
    print("S(R) = exp(-(ζ*R)/2) * (1 + (1/2)*(ζ*R) + (1/12)*(ζ*R)**2 + (1/240)*(ζ*R)**4)")


# Execute the function to print the result
print_overlap_integral_2s_2s()
<<<S(R) = exp(-(ζ*R)/2) * (1 + (1/2)*(ζ*R) + (1/12)*(ζ*R)**2 + (1/240)*(ζ*R)**4)>>>