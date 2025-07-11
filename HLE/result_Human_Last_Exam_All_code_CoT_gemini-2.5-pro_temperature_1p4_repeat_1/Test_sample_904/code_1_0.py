import sympy

def find_coefficients():
    """
    Symbolically determines and prints the coefficients A(r) and B(r) for the
    governing linear equation of the fluid interface shape xi(r).
    
    The final equation has the form:
    A(r) * d^2(xi)/dr^2 + B(r) * d(xi)/dr + C(r, xi) = 0
    """
    
    # Define symbolic variables used in the problem
    r = sympy.Symbol('r')
    gamma = sympy.Symbol('gamma')  # Represents the surface tension
    xi = sympy.Function('xi')(r)   # Represents the interface displacement xi as a function of r

    # Define the first and second derivatives of xi(r)
    xi_prime = xi.diff(r)
    xi_double_prime = xi.diff(r, 2)
    
    # The linearized governing equation is derived from the Young-Laplace equation:
    # gamma * (d^2(xi)/dr^2 + (1/r) * d(xi)/dr) - P_el(r, xi) = 0
    # where P_el is the electrostatic pressure.
    
    # We construct the part of the equation involving the derivatives of xi
    equation_lhs = gamma * xi_double_prime + (gamma / r) * xi_prime
    
    # Use sympy's 'coeff' method to extract the coefficients of the derivatives
    # A(r) is the coefficient of the second derivative
    A_r = equation_lhs.coeff(xi_double_prime)
    
    # B(r) is the coefficient of the first derivative
    B_r = equation_lhs.coeff(xi_prime)
    
    # Print the final results for the coefficients
    print("The governing linear equation for the interfacial shape xi(r) is:")
    print("gamma * d^2(xi)/dr^2 + (gamma/r) * d(xi)/dr - P_el(r, xi) = 0")
    print("-" * 50)
    print("Comparing this to the form A(r)*xi'' + B(r)*xi' + C(r, xi) = 0, we find:")
    print("\nThe coefficient A(r) is:")
    print(f"A(r) = {A_r}")
    print("\nThe coefficient B(r) is:")
    print(f"B(r) = {B_r}")

if __name__ == "__main__":
    find_coefficients()
