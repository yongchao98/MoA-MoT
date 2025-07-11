import math

def print_overlap_integral_formula():
    """
    Prints the analytical expression for the overlap integral of two 2s
    hydrogenic orbitals.
    """
    
    # The derived analytical expression for the overlap integral S.
    # The variable 'rho' in the formula represents the product of the effective
    # nuclear charge (zeta) and the internuclear distance (R), i.e., rho = zeta * R.

    print("The analytical expression for the overlap integral S for two 2s hydrogenic orbitals is:")
    print("\nLet rho = zeta * R, where 'zeta' is the effective nuclear charge and 'R' is the internuclear distance.")
    print("The overlap integral S(rho) is then given by:")

    # Print the equation with its terms clearly laid out
    print("\nS(rho) = exp(-rho / 2) * (1 + (rho / 2) + (rho**2 / 12) + (rho**4 / 240))")
    
    print("\nLet's break down the coefficients in the polynomial term:")
    print(f"Term with rho^0: 1 = {1}")
    print(f"Coefficient for rho^1: 1/2 = {1/2}")
    print(f"Coefficient for rho^2: 1/12 = {1/12:.6f}")
    # Note: There is no rho^3 term in the final expression.
    print(f"Coefficient for rho^3: 0 = {0}")
    print(f"Coefficient for rho^4: 1/240 = {1/240:.6f}")


if __name__ == '__main__':
    print_overlap_integral_formula()
