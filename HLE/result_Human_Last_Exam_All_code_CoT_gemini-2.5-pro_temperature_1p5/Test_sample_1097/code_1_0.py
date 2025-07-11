import sympy as sp

def solve():
    """
    This script defines and justifies the set M for the Banach Fixed-Point Theorem
    to solve the given Boundary Value Problem.
    """
    
    # Define symbols
    x, s = sp.symbols('x s', real=True, positive=True)

    print("Step 1: Reformulate BVP as a fixed-point problem u = T(u)")
    print("---------------------------------------------------------")
    print("The BVP is u''(x) = exp(u(x)), with u(0)=u(1)=0.")
    print("The operator T is (Tu)(x) = integral from 0 to 1 of G(x, s) * exp(u(s)) ds.")
    print("where G(x, s) is the Green's function for u'' with zero boundary conditions on [0, 1].\n")

    print("Step 2: Define the set M")
    print("---------------------------------------------------------")
    print("Analysis shows any solution u(x) must be convex (since u''>0) and satisfy u(0)=u(1)=0.")
    print("This implies that any solution must be non-positive, i.e., u(x) <= 0.")
    print("Therefore, the correct choice for M is the set of non-positive continuous functions.")
    
    m_definition = "M = {u in C([0, 1]) | u(x) <= 0 for all x in [0, 1]}"
    print(f"The proposed set is: {m_definition}\n")

    print("Step 3: Verify T is a contraction on M")
    print("---------------------------------------------")
    print("To prove T is a contraction, we need to show ||Tu - Tv|| <= k * ||u - v|| for k < 1.")
    print("The contraction constant k is given by: k = max_x integral_0^1 |G(x,s)| ds.")
    
    # G(x,s) is non-positive, so |G(x,s)| = -G(x,s).
    # G(x,s) = x(s-1) for s < x and s(x-1) for s > x.
    # Note: Sympy's integrate has issues with piecewise functions defined with inequalities.
    # We split the integral manually.
    # We must use different symbols for integration limits vs. variables inside
    t = sp.Symbol('t')
    integral_part1 = sp.integrate(-x * (t - 1), (t, x, 1))
    integral_part2 = sp.integrate(-t * (x - 1), (t, 0, x))
    integral_G = sp.simplify(integral_part1 + integral_part2)

    print(f"The integral of |G(x,s)| with respect to s from 0 to 1 is: {integral_G}")

    # To find the maximum of x/2 - x**2/2 on [0,1], we find the critical point.
    derivative = sp.diff(integral_G, x)
    # The critical point is where the derivative is 0, which is at x = 1/2.
    critical_point = sp.solve(derivative, x)[0]
    
    # The maximum value is at this critical point.
    k_value = integral_G.subs(x, critical_point)

    print(f"The maximum value of this expression occurs at x = {critical_point}, giving k = {k_value}.")
    print(f"Since k = {k_value} < 1, T is a contraction on M.\n")
    
    print("Conclusion")
    print("----------")
    print("Since M is a complete metric space and T is a contraction from M to M, the Banach Fixed-Point Theorem applies.")
    print("The final answer is the definition of the set M.\n")

    print(f"The set M you should define is:\n{m_definition}")

solve()