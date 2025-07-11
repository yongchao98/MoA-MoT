import sympy as sp

def verify_separatrix():
    """
    This function verifies if the curve d = -u**2 is a separatrix for the given system of ODEs.
    A separatrix is an invariant manifold, meaning a trajectory that starts on it stays on it.
    We verify this by substituting d = -u**2 into the first ODE and checking if the equation holds.
    """
    # Define time t as a symbol and u, d as functions of t
    t = sp.Symbol('t')
    u = sp.Function('u')(t)
    d = sp.Function('d')(t)

    # Define the system of differential equations
    # d'(t) = 2*d(t)**2 + (-3*u(t) + 5*u(t)**2)*d(t) - u(t)**3*(1-u(t))
    d_prime_expr = 2*d**2 + (-3*u + 5*u**2)*d - u**3*(1-u)
    # u'(t) = u(t)**2*(u(t)-1)
    u_prime_expr = u**2 * (u - 1)

    # Define the candidate separatrix equation: d(t) = -u(t)**2
    separatrix_candidate = -u**2

    # Calculate the left-hand side (LHS) of the first ODE: d'(t)
    # We differentiate our candidate d with respect to t
    lhs = sp.diff(separatrix_candidate, t)
    # The result contains u'(t), which we substitute with its expression from the second ODE
    lhs = lhs.subs(sp.diff(u, t), u_prime_expr)

    # Calculate the right-hand side (RHS) of the first ODE
    # We substitute our candidate for d into the expression for d'(t)
    rhs = d_prime_expr.subs(d, separatrix_candidate)

    # Check if LHS and RHS are equal by simplifying their difference.
    # If the difference simplifies to 0, the candidate is a valid invariant manifold.
    if sp.simplify(lhs - rhs) == 0:
        print("Verification successful.")
        print("The equation d = -u**2 is a separatrix of the system.")
        
        # As requested, output the numbers in the final equation.
        # The equation is of the form d = a * u**p
        a = -1
        p = 2
        
        print("\nThe equation for the separatrix is d = a * u**p")
        print(f"The coefficient a is: {a}")
        print(f"The power p is: {p}")
        print("\nFinal Equation:")
        print(f"d = ({a}) * u**({p})")

    else:
        print("Verification failed.")
        print("The curve d = -u**2 is not a separatrix.")

# Run the verification
verify_separatrix()