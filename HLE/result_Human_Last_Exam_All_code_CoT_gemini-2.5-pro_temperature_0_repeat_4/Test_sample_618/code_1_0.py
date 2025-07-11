def print_transformation_on_x():
    """
    This function prints the general representation for the Lie group transformation
    on the spatial coordinate 'x' for the heat equation with a logarithmic source term:
    u_t = u_{xx} + (k1*ln(u) + k2)*u.
    """
    
    print("The infinitesimal generator for the transformation on the spatial coordinate x is denoted by xi.")
    print("Solving the determining equations yields:")
    print("xi = C4 - (2 * C5 / k1) * exp(k1 * t)")
    print("where C4 and C5 are arbitrary constants, and k1 is the parameter from the PDE.\n")

    print("The finite transformation x' = g(t, x; lambda) is found by integrating the ODE d(x')/d(lambda) = xi(t', x').")
    print("The result depends on the time transformation infinitesimal, tau = C3, leading to two cases.\n")
    
    print("="*70)
    print("Case 1: The time transformation constant C3 is non-zero (C3 != 0)")
    print("="*70)
    print("The finite transformation on x is:")
    print("x' = x + lambda*C4 - (2*C5 / (k1**2 * C3)) * exp(k1*t) * (exp(k1*C3*lambda) - 1)")
    print("\nWhere:")
    print("  x'      : The transformed spatial coordinate")
    print("  x       : The original spatial coordinate")
    print("  t       : The time variable")
    print("  lambda  : The continuous group parameter")
    print("  k1      : The constant from the PDE's source term")
    print("  C3, C4, C5: Arbitrary constants defining a specific transformation\n")

    print("="*70)
    print("Case 2: The time transformation constant C3 is zero (C3 = 0)")
    print("="*70)
    print("The finite transformation on x simplifies to:")
    print("x' = x + lambda * (C4 - (2*C5 / k1) * exp(k1*t))")
    print("\nWhere the variables and constants have the same meaning as in Case 1.\n")

if __name__ == '__main__':
    print_transformation_on_x()