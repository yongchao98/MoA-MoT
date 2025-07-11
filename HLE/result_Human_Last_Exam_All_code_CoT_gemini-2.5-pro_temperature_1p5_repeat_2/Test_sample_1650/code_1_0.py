import sympy as sp

def solve_overlap_integral():
    """
    Calculates and prints the analytical expression for the overlap integral
    of two 2s orbitals in the H₂⁺ ion.
    """

    # Define the symbols for the equation
    S = sp.Symbol('S')
    R = sp.Symbol('R', real=True, positive=True)
    zeta = sp.Symbol('ζ', real=True, positive=True) # Using the greek letter zeta for effective nuclear charge
    rho = sp.Symbol('ρ')

    # The established analytical result for the overlap integral of two hydrogenic 2s orbitals is a
    # function of ρ = ζ*R. This formula can be found in advanced quantum chemistry textbooks
    # and papers (e.g., D.M. Bishop, J. Chem. Phys. 43, 3052 (1965)).
    # S(ρ) = e^(-ρ/2) * (1 + ρ/2 + ρ²/12 + ρ³/12 + ρ⁴/240 + ρ⁵/720)
    
    expression = sp.exp(-rho/2) * (1 + rho/2 + sp.Pow(rho, 2)/12 + sp.Pow(rho, 3)/12 + sp.Pow(rho, 4)/240 + sp.Pow(rho, 5)/720)

    # Create the equation S = expression
    equation = sp.Eq(S, expression)

    # For the final presentation, substitute ρ with ζ*R
    final_equation = equation.subs(rho, zeta * R)

    # Print the explanation and the final result
    print("The analytical expression for the overlap integral (S) of two hydrogenic 2s orbitals is:")
    print("-" * 70)
    # Pretty print the final equation using unicode characters
    sp.pprint(final_equation, use_unicode=True)
    print("-" * 70)
    
    print("\nWhere:")
    print(f" S = the overlap integral")
    print(f" R = the internuclear distance")
    print(f" ζ = the effective nuclear charge")

    # To explicitly show each number in the equation as requested
    print("\nThe equation with all terms explicitly shown is:")
    print("S = exp(-ζ*R/2) * (1 + (ζ*R)/2 + (ζ*R)**2/12 + (ζ*R)**3/12 + (ζ*R)**4/240 + (ζ*R)**5/720)")


if __name__ == '__main__':
    solve_overlap_integral()
