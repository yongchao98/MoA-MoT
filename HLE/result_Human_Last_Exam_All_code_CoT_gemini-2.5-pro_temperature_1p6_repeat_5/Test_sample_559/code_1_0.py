def find_separatrix():
    """
    This function presents the solution for the separatrix of the given system of ODEs.

    The analytical derivation shows that the separatrix is an invariant manifold
    that separates the basin of attraction of the origin from a region of solutions
    that blow up. This curve is found to have the equation:
    d = u - u^2

    This equation can be rewritten in the standard form c1*d + c2*u + c3*u^2 = 0.
    """
    
    # Coefficients for the separatrix equation: d - u + u^2 = 0
    c1 = 1  # Coefficient for d
    c2 = -1 # Coefficient for u
    c3 = 1  # Coefficient for u^2
    
    print("The separatrix for the given system of differential equations has been found analytically.")
    print("Its equation is:")
    print("d = u - u^2")
    print("\nThis can be written in the form c1*d + c2*u + c3*u^2 = 0.")
    print(f"The coefficients are c1 = {c1}, c2 = {c2}, c3 = {c3}.")
    
    print("\nThe full equation with each numerical term is:")
    print(f"{c1} d(t) + ({c2}) u(t) + {c3} u(t)^2 = 0")

find_separatrix()