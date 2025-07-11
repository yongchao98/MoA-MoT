import sympy as sp

def define_set_m_for_bvp():
    """
    This function performs the symbolic calculations needed to define the set M
    for the Banach Fixed-Point Theorem applied to the given BVP.
    """
    # Define symbols for our variables
    x, s = sp.symbols('x s')

    print("Step 1: Finding the properties of the integral operator.")
    # The BVP u''(x) = exp(u(x)), u(0)=u(1)=0 is equivalent to the integral equation
    # u(x) = integral(G(x,s) * exp(u(s)), s=0..1), where G(x,s) is the Green's function.
    # The Green's function is G(x,s) = s(x-1) for s<=x and x(s-1) for s>=x.
    # Note that G(x,s) <= 0 for x, s in [0,1].
    # To check the contraction and self-mapping properties, we need to evaluate the
    # maximum value of the integral of |G(x,s)| with respect to s.
    # Let h(x) = integral(|G(x,s)|, s=0..1).
    # Since G(x,s) is negative, |G(x,s)| = -G(x,s).
    # |G(x,s)| = s(1-x) for s<=x and x(1-s) for s>=x.

    integral_part1 = sp.integrate(s * (1 - x), (s, 0, x))
    integral_part2 = sp.integrate(x * (1 - s), (s, x, 1))
    h_x = sp.simplify(integral_part1 + integral_part2)

    print(f"The integral h(x) = integral(|G(x,s)|, s=0..1) is: {h_x}")

    # To find the maximum of h(x) on the interval [0,1], we find its derivative.
    h_prime = sp.diff(h_x, x)
    
    # Solve for critical points where the derivative is zero.
    critical_points = sp.solve(h_prime, x)
    critical_point = critical_points[0]

    # The maximum value of the integral occurs at this critical point.
    max_integral_val = h_x.subs(x, critical_point)

    print(f"The maximum value of h(x) on [0, 1] is K = {max_integral_val}")

    print("\nStep 2: Defining the set M based on the calculations.")
    # The contraction constant for the operator T is K * Lipschitz_const(exp(u)).
    # Since any solution u(x) <= 0, we can define M as a set of non-positive functions.
    # For u <= 0, the Lipschitz constant for exp(u) is bounded by exp(0) = 1.
    # So, the contraction constant is K * 1 = 1/8, which is < 1.
    
    # For the self-mapping property T(M) subset M, we need ||Tu|| <= R for u in M.
    # ||Tu|| <= max_integral_val * sup|exp(u)| = 1/8 * 1 = 1/8.
    # So, we need to choose the bound R of our set M such that R >= 1/8.
    # The simplest choice is R = 1/8.

    R = max_integral_val
    u_lower_bound = -R
    u_upper_bound = 0
    x_lower_bound = 0
    x_upper_bound = 1

    print("\n--- Final Definition of the Set M ---")
    print("The set M is the collection of continuous functions u(x) defined on the interval "
          f"[{x_lower_bound}, {x_upper_bound}] that satisfy:")
    print(f"1. Boundary conditions: u({x_lower_bound}) = {u_upper_bound} and u({x_upper_bound}) = {u_upper_bound}.")
    print(f"2. A uniform bound: {u_lower_bound} <= u(x) <= {u_upper_bound} for all x in the interval.")
    
define_set_m_for_bvp()