def find_separatrix_equation():
    """
    This function describes the separatrix of the given system of ODEs.
    The separatrix is found to be the curve d = -u**2, or d + u**2 = 0.
    This script will print the equation in a fully explicit form.
    """

    # The equation of the separatrix is d + u**2 = 0.
    # We can write this in a general form c1*d**p1 + c2*u**p2 = c3
    # where c are coefficients and p are powers.
    d_coeff = 1
    d_power = 1
    u_coeff = 1
    u_power = 2
    rhs_constant = 0

    # Print the equation with all numbers explicitly shown.
    print("The equation for the separatrix is:")
    print(f"{d_coeff} * d**{d_power} + {u_coeff} * u**{u_power} = {rhs_constant}")

find_separatrix_equation()