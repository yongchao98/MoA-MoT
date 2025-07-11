import sympy as sp

def display_overlap_integral_expression():
    """
    Derives and displays the analytical expression for the overlap integral
    between two 2s orbitals in the H2+ ion.
    """
    # Define symbolic variables
    R, zeta = sp.symbols('R zeta', real=True, positive=True)
    rho = sp.Symbol('rho', real=True, positive=True)

    # Final derived expression in terms of rho
    S_rho = sp.exp(-rho) * (1 + rho + sp.Rational(1, 3) * rho**2 + sp.Rational(1, 15) * rho**4)

    # Substitute rho = zeta*R/2 to get the expression in terms of R and zeta
    S_R_zeta = S_rho.subs(rho, zeta * R / 2)

    # --- Output Section ---
    print("This script provides the analytical expression for the overlap integral S.\n")
    print("Let R be the internuclear distance and zeta be the effective nuclear charge.")
    print("A convenient parameter, rho, is defined as rho = (zeta * R) / 2.\n")
    
    # Print expression in terms of rho
    print("The final expression for the overlap integral in terms of rho is:")
    print("S(rho) = exp(-rho) * (", end="")
    print("1", "+", "rho", "+", "(", "1", "/", "3", ")*rho**2", "+", "(", "1", "/", "15", ")*rho**4", ")", sep=' ')
    print("-" * 50)

    # Print expression in terms of R and zeta
    print("The final expression in terms of R and zeta is:")
    
    # Extract coefficients for printing
    poly_R_zeta = sp.expand(S_R_zeta.args[1])
    c0 = poly_R_zeta.coeff(zeta, 0).coeff(R, 0)
    c1 = poly_R_zeta.coeff(zeta**1 * R**1)
    c2 = poly_R_zeta.coeff(zeta**2 * R**2)
    c3 = poly_R_zeta.coeff(zeta**3 * R**3) # This will be 0
    c4 = poly_R_zeta.coeff(zeta**4 * R**4)

    print(f"S(R, zeta) = exp(-zeta*R / 2) * (", end="")
    print(f"{c0}", "+", f"({c1})*zeta*R", end="")
    if c2 != 0:
        print(" +", f"({c2})*zeta**2*R**2", end="")
    if c3 != 0: # This term is zero in the final result, but included for robustness
        print(" +", f"({c3})*zeta**3*R**3", end="")
    if c4 != 0:
        print(" +", f"({c4})*zeta**4*R**4", end="")
    print(" )")

if __name__ == '__main__':
    display_overlap_integral_expression()