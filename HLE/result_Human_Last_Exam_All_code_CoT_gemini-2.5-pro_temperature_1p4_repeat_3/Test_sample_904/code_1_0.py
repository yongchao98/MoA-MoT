import sympy

def solve_fluid_equation_coeffs():
    """
    This function derives and prints the coefficients A(r) and B(r) for the governing
    linear equation of the fluid interface shape.
    """
    # Define the symbols used in the equation
    r, gamma = sympy.symbols('r gamma')

    # From the derivation based on the linearized Young-Laplace equation in cylindrical coordinates,
    # the governing equation is:
    # gamma * d^2(xi)/dr^2 + (gamma/r) * d(xi)/dr - P_elec(r) = 0
    #
    # Comparing this to the standard form:
    # A(r) * d^2(xi)/dr^2 + B(r) * d(xi)/dr + C(r, xi) = 0

    # Identify the coefficients A(r) and B(r)
    A_r = gamma
    B_r = gamma / r

    # Print the results
    print("The governing linear equation for the interfacial shape is of the form:")
    print("A(r) * d^2(xi)/dr^2 + B(r) * d(xi)/dr + C(r, xi) = 0")
    print("\nBased on the derivation from the linearized Young-Laplace equation:")
    print(f"The coefficient A(r) is: {A_r}")
    print(f"The coefficient B(r) is: {B_r}")

solve_fluid_equation_coeffs()