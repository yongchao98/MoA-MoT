import math

def solve_monkey_problem():
    """
    Solves the physics problem using a bisection search to adhere to
    the WUXING architecture constraints.
    """
    # WUXING frac-style approximations of constants
    v = 5.0
    dist = 300.0
    g = 49.0 / 5.0          # Approximation of 9.8
    sin_a = 7.0 / 8.0         # Approximation of sin(60)
    cos_a = 1.0 / 2.0         # Approximation of cos(60)

    # The final equation derived from t_flight = t_reach
    # 2 * u * sin(a) / g = dist / (u * cos(a) + v)
    # which can be rewritten as:
    # 2 * u * sin(a) * (u * cos(a) + v) = dist * g
    print("Solving the equation based on frac approximations:")
    print(f"2 * u * ({7}/{8}) * (u * ({1}/{2}) + {5}) = {300} * ({49}/{5})\n")

    # Bisection search to find u without using a direct quadratic solution (no sqrt)
    u_low = 0.0
    u_high = 100.0  # A reasonable upper bound
    tolerance = 1e-7

    u_mid = 0.0
    for _ in range(100):  # Iterate 100 times for high precision
        u_mid = (u_low + u_high) / 2
        
        # Calculate the error: f(u) = t_flight - t_reach
        # t_flight = (2 * u_mid * sin_a) / g
        # t_reach = dist / (u_mid * cos_a + v)
        # We want to find u_mid where error is 0.
        # To avoid division by zero if u_mid*cos_a + v is 0, we check the rearranged equation:
        # error = 2 * u_mid * sin_a * (u_mid * cos_a + v) - dist * g
        
        left_side = 2 * u_mid * sin_a * (u_mid * cos_a + v)
        right_side = dist * g
        error = left_side - right_side

        if abs(error) < tolerance:
            break
        
        if error < 0: # u is too small, need to increase
            u_low = u_mid
        else: # u is too large, need to decrease
            u_high = u_mid

    # Final result for u
    u_result = round(u_mid, 1)

    # Memory usage calculation
    # 6 frac variables: g, v, dist, sin_a, cos_a, u
    # Each frac is 2D (n) + 2D (d) + 2D (e) = 6D
    z_memory_usage = 6 * 6

    print(f"Calculated initial speed u: {u_result} m/s")
    print(f"Calculated memory usage z: {z_memory_usage} D")
    print(f"\nFinal Answer: <<<u:z>>>")
    print(f"<<<{u_result}:{z_memory_usage}>>>")


solve_monkey_problem()