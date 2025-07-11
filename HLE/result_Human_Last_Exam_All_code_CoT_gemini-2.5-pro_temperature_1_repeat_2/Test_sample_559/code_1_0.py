def find_separatrix():
    """
    This function presents the equation of the separatrix for the given system of ODEs.

    The analysis shows that the separatrix is the unstable manifold of the saddle
    point at (u, d) = (1, -1). This manifold is described by the algebraic
    equation d = -u^2.
    """

    # The equation of the separatrix is d = -u^2.
    # To meet the requirement of showing each number in the equation,
    # we can write it in the form: a*d + b*u^c = k
    # This gives: 1*d + 1*u^2 = 0

    coeff_d = 1
    coeff_u_term = 1
    power_u = 2
    constant = 0

    print("The equation of the separatrix is d = -u^2")
    print("\nTo explicitly output each number in the final equation:")
    print(f"Equation form: a*d + b*u^c = k")
    print(f"The equation is: {coeff_d}*d + {coeff_u_term}*u^{power_u} = {constant}")

find_separatrix()