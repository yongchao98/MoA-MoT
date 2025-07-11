import sympy

def find_coefficients():
    """
    This function determines and prints the coefficients A(r) and B(r) for the
    governing linear equation of the fluid interface shape xi(r), based on the
    linearized Young-Laplace equation.
    """
    # Define symbols for the variables. 'gamma' for surface tension and 'r' for radial position.
    gamma, r = sympy.symbols('gamma r')

    # From the derivation, the linearized pressure balance equation is:
    # gamma * xi''(r) + (gamma/r) * xi'(r) + C(r, xi) = 0
    #
    # We compare this to the general form given in the problem:
    # A(r) * xi''(r) + B(r) * xi'(r) + C(r, xi) = 0

    # The coefficient of the second derivative term xi''(r)
    A_r = gamma

    # The coefficient of the first derivative term xi'(r)
    B_r = gamma / r

    print("The governing linear equation for the interfacial shape xi(r) is of the form:")
    print("A(r) * d^2(xi)/dr^2 + B(r) * d(xi)/dr + C(r, xi) = 0\n")
    print("Based on the linearization of the Young-Laplace equation in cylindrical coordinates, the coefficients are:\n")

    # Output each component of the final equation
    print(f"A(r) = {A_r}")
    print(f"B(r) = {B_r}")

# Execute the function to display the results
find_coefficients()