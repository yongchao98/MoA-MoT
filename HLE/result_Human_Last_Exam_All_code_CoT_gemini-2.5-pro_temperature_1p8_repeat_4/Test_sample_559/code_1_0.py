import sympy

def solve():
    """
    This function verifies the separatrix for the given system of ODEs and prints the result.
    """
    # Define symbols for time t, and functions u(t) and d(t)
    t = sympy.Symbol('t')
    u = sympy.Function('u')
    d = sympy.Function('d')

    # Define the system of differential equations
    # d'(t) = 2d^2 + (5u^2-3u)d - u^3(1-u)
    # u'(t) = u^2(u-1)
    ode_d_rhs = 2*d(t)**2 + (5*u(t)**2 - 3*u(t))*d(t) - u(t)**3 * (1-u(t))
    ode_u_rhs = u(t)**2 * (u(t)-1)

    # The equation for u'(t)
    u_prime_eq = sympy.Eq(sympy.Derivative(u(t), t), ode_u_rhs)

    # We propose the separatrix is given by the curve d = -u^2.
    # Let's verify if d(t) = -u(t)**2 is a valid trajectory.
    proposed_d = -u(t)**2

    # Calculate the Left-Hand Side (LHS) of the first ODE, which is d'(t).
    # We use the chain rule: d/dt = (d/du)*(du/dt)
    # So, d/dt(-u(t)**2) = -2*u(t)*u'(t)
    lhs = sympy.Derivative(proposed_d, t).doit()
    # Substitute u'(t) from the second ODE
    lhs = lhs.subs(sympy.Derivative(u(t), t), u_prime_eq.rhs)

    # Calculate the Right-Hand Side (RHS) of the first ODE by substituting d = -u^2.
    rhs = ode_d_rhs.subs(d(t), proposed_d)

    # Check if LHS and RHS are equal. sympy.simplify helps to confirm this.
    if sympy.simplify(lhs - rhs) == 0:
        # The separatrix is d = -u^2. We present it in the form d = a * u^b.
        a = -1
        b = 2

        print("The separatrix is given by the equation d = -u^2.")
        print("This has been symbolically verified.")
        print("The equation is of the form: d = a * u^b")
        # Outputting each number in the final equation
        print(f"The coefficient 'a' is: {a}")
        print(f"The exponent 'b' is: {b}")
        print(f"Final Equation: d = {a}u^{b}")
    else:
        print("Verification failed. The proposed solution is not correct.")

solve()