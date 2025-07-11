import sympy

def find_separatrix():
    """
    This function verifies the proposed separatrix equation and prints its components.
    """
    # Define u as a symbol
    u = sympy.Symbol('u')

    # The equation of the separatrix we found through analysis is d = -u^2.
    # Let's define it.
    d = -u**2

    # Now, let's define the original system's trajectory equation dd/du.
    # dd/du = (d'(t)) / (u'(t))
    dd_dt = 2*d**2 + (5*u**2 - 3*u)*d - u**3*(1-u)
    du_dt = u**2*(u - 1)
    
    # Avoid division by zero at the equilibria
    if u_dt == 0:
        # This check is symbolic, so we look for simplification.
        # This part of the code is for verification.
        pass
    else:
        slope_from_system = dd_dt / du_dt

    # The slope can also be found by differentiating our solution d(u) with respect to u.
    slope_from_solution = sympy.diff(d, u)

    # To verify our solution, we check if the two slopes are identical.
    # We do this by simplifying their difference; it should be zero.
    difference = sympy.simplify(slope_from_system - slope_from_solution)

    # Print the verification result.
    if difference == 0:
        print("Verification successful: The equation d = -u**2 is an exact solution for a trajectory of the system.")
    else:
        print("Verification failed. The proposed solution is not correct.")
        print(f"The simplified difference is: {difference}")

    # The equation of the separatrix is d = -1 * u**2.
    # As requested, we will output each number in the final equation.
    # The numbers are the coefficient (-1) and the exponent (2).
    c = -1
    p = 2
    print("\nThe separatrix is of the form d = c * u**p")
    print(f"The coefficient c is: {c}")
    print(f"The exponent p is: {p}")

find_separatrix()