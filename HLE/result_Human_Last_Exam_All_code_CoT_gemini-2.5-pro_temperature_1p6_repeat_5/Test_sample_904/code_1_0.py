import sympy

def solve_fluid_equation():
    """
    This function defines the physical parameters as symbolic variables and
    prints the derived coefficients A(r) and B(r) for the governing
    differential equation of the fluid interface.
    """
    # Define the symbols for surface tension (gamma) and radial position (r)
    gamma, r = sympy.symbols('gamma r', positive=True, real=True)
    xi = sympy.Function('xi')(r)

    # The governing differential equation for the interface xi(r) is derived
    # from the balance of surface tension and pressure forces.
    # The linearized form of this equation is:
    # (gamma * r) * xi''(r) + gamma * xi'(r) - r * Delta_P(r, xi) = 0
    #
    # Comparing this to the general form A(r)*xi''(r) + B(r)*xi'(r) + C(r,xi) = 0,
    # we identify the coefficients A(r) and B(r).

    # Coefficient of the second derivative term
    A_r = gamma * r

    # Coefficient of the first derivative term
    B_r = gamma

    print("The governing linear differential equation for the interfacial shape xi(r) is of the form:")
    print("A(r) * d^2(xi)/dr^2 + B(r) * d(xi)/dr + C(r, xi) = 0")
    print("\nBased on the derivation from the force balance on an interface element:\n")
    
    # We use sympy's pretty print for a clear mathematical output
    print("Coefficient A(r):")
    sympy.pprint(A_r)
    print("\nCoefficient B(r):")
    sympy.pprint(B_r)

if __name__ == "__main__":
    solve_fluid_equation()
