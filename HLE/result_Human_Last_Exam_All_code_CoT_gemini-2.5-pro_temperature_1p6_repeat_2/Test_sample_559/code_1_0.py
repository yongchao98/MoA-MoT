def find_separatrix_equations():
    """
    This script finds and prints the equations for the separatrices of the given system.

    The system of differential equations is:
    d'(t) = 2*d(t)^2 + (-3*u(t) + 5*u(t)^2)*d(t) - u(t)^3*(1-u(t))
    u'(t) = u(t)^2*(u(t)-1)

    The method involves finding the equilibrium points and identifying the saddle point.
    The separatrices are the stable and unstable manifolds of this saddle point.
    Through mathematical analysis, the saddle point is located at (u, d) = (1, -1).
    The equations for its stable and unstable manifolds are then determined.
    """
    
    print("The system has a saddle point at (u, d) = (1, -1).")
    print("The separatrices are the stable and unstable manifolds of this point.")
    print("-" * 30)

    # The first separatrix is the stable manifold of the saddle point.
    # Mathematical analysis shows this is a vertical line.
    print("One separatrix is the line given by the equation:")
    u_coeff_1 = 1
    const_1 = 1
    print(f"The equation is: {u_coeff_1} * u = {const_1}")
    print("\nBreaking down the numbers in this equation:")
    print(f"The coefficient of the variable u is {u_coeff_1}.")
    print(f"The constant value on the right-hand side is {const_1}.")
    
    print("\n" + "-" * 30 + "\n")
    
    # The second separatrix is the unstable manifold of the saddle point.
    # Mathematical analysis shows this is a parabolic curve.
    print("The other separatrix is the curve given by the equation:")
    d_coeff_2 = 1
    sign_2 = "-"
    u_coeff_2 = 1
    u_power_2 = 2
    print(f"The equation is: {d_coeff_2} * d = {sign_2}{u_coeff_2} * u^{u_power_2}")
    print("\nBreaking down the numbers in this equation:")
    print(f"The coefficient of the variable d is {d_coeff_2}.")
    print(f"The sign of the term on the right is negative, corresponding to a multiplicative factor of -{u_coeff_2}.")
    print(f"The coefficient of the u-term is {u_coeff_2}.")
    print(f"The power of u is {u_power_2}.")

# Run the function to display the results.
find_separatrix_equations()