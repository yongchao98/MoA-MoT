import sympy

def solve_separatrix():
    """
    This function analyzes the given system of ODEs to find the separatrix.
    The analysis is done analytically, and the code prints the final result.
    """

    # The system of differential equations is:
    # d'(t) = 2*d^2(t) + (-3*u(t) + 5*u^2(t))*d(t) - u^3(t)*(1-u(t))
    # u'(t) = (u(t) - 1)*u^2(t)

    # Step 1 & 2: Find and classify equilibrium points.
    # From u'(t) = 0, we get u=0 or u=1.
    # For u=0, d'(t)=0 gives 2*d^2 = 0, so d=0. Point: (0, 0).
    # For u=1, d'(t)=0 gives 2*d^2 + 2*d = 0, so d=0 or d=-1. Points: (1, 0), (1, -1).
    # Analysis of the Jacobian matrix shows that (1, -1) is a saddle point.

    # Step 3 & 4: Find an invariant manifold for the saddle point (1, -1).
    # A separatrix is often the stable or unstable manifold of a saddle point.
    # We test if a curve of the form d = f(u) is an invariant manifold.
    # By testing simple polynomial forms, we can discover that the curve d = -u^2 is an invariant manifold.
    # Let's verify this. Let H(u, d) = d + u^2. An invariant manifold H=0 exists if dH/dt = 0 whenever H=0.
    # dH/dt = d(d)/dt + d(u^2)/dt = d' + 2*u*u'
    # Substituting the expressions for d' and u':
    # dH/dt = (2*d^2 + (5*u^2-3*u)*d - u^3*(1-u)) + 2*u*(u^2*(u-1))
    # dH/dt = 2*d^2 + (5*u^2-3*u)*d + 3*u^4 - 3*u^3
    # Now, we evaluate dH/dt on the curve d = -u^2:
    # dH/dt = 2*(-u^2)^2 + (5*u^2-3*u)*(-u^2) + 3*u^4 - 3*u^3
    #        = 2*u^4 - 5*u^4 + 3*u^3 + 3*u^4 - 3*u^3 = 0.
    # This confirms that d = -u^2 is an invariant curve.

    # Step 5: This curve is the separatrix.
    # The curve d = -u^2 passes through the saddle point (1, -1) because -1 = -(1)^2.
    # It can be shown this is the unstable manifold of the saddle point.
    # Therefore, it is a separatrix of the system.

    print("The separatrix of the system is given by the equation of a parabola.")
    
    # Define symbols for pretty printing
    d, u = sympy.symbols('d u')
    
    # Define the equation
    # Coefficients for the equation in the form d = a*u**2 + b*u + c
    a = -1
    b = 0
    c = 0
    
    # The equation is d = -u**2
    equation = sympy.Eq(d, a*u**2 + b*u + c)

    print(f"The equation is: {sympy.pretty(equation)}")
    
    # Printing each number in the equation explicitly.
    print("\nBreaking down the equation d = -u^2:")
    print("The equation can be written in the form: 1*d = a*u^b + ...")
    print(f"The coefficient of d is: 1")
    print(f"The coefficient 'a' of u^2 is: {a}")
    print(f"The exponent 'b' of u is: 2")
    
solve_separatrix()
