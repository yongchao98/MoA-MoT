import math

def solve_projectile_speed():
    """
    This function calculates the initial speed 'u' for the rock to hit the lion
    and the memory usage 'z' based on the Wuxing architecture constraints.
    """
    # 1. Define physical constants and Wuxing-specific approximations.
    d_lion = 300.0  # Initial distance in meters
    v_lion = 5.0    # Lion's speed in m/s
    angle_deg = 60.0 # Throw angle in degrees

    # In the Wuxing architecture, g=10 is assumed for simpler decimal arithmetic.
    g = 10.0

    # The Wuxing architecture lacks sin() and sqrt(), so we use rational approximations
    # that fit into the 'frac' type's 'char' components.
    # sqrt(3) is approximated as 43/25 = 1.72.
    sqrt3_approx = 43.0 / 25.0
    sin_a = sqrt3_approx / 2.0
    cos_a = 1.0 / 2.0

    # 2. Derive coefficients for the quadratic equation: A*u^2 + B*u + C = 0
    # Equation form: (sin(2a))*u^2 + (2*v*sin(a))*u - d*g = 0
    A = 2 * sin_a * cos_a
    B = 2 * v_lion * sin_a
    C = -d_lion * g

    # Output the final equation with numerical values as required.
    print(f"Derived equation for u: ({A:.2f}) * u^2 + ({B:.2f}) * u + ({C:.0f}) = 0")

    # 3. Solve the equation using an iterative method (Newton-Raphson)
    # as the sqrt() function is unavailable on the Wuxing machine.
    # f(u) = A*u^2 + B*u + C
    # f'(u) = 2*A*u + B
    # u_new = u_old - f(u) / f'(u)

    u = 50.0  # Initial guess for the speed u.

    # Iterate a few times to find a precise solution.
    for _ in range(5):
        f_u = A * u**2 + B * u + C
        f_prime_u = 2 * A * u + B
        if f_prime_u == 0:
            break
        u = u - f_u / f_prime_u

    # Round the final speed to one decimal place.
    u_final = round(u, 1)

    # 4. Calculate memory usage 'z'.
    # We need 6 primary variables: d_lion, v_lion, g, sin_a, cos_a, u.
    # Each 'frac' variable occupies 3 * 2D = 6D.
    # Total memory z = 6 variables * 6D/variable = 36D.
    z = 36

    # 5. Print the final answer in the specified format 'u:z'.
    print(f"\nCalculated initial speed u: {u_final} m/s")
    print(f"Calculated memory usage z: {z} D")
    print(f"Final Answer (u:z) -> {u_final}:{z}")

# Execute the solution function.
solve_projectile_speed()

<<<53.9:36>>>