def solve_monkey_and_lion():
    """
    This function solves the physics problem based on the Wuxing architecture constraints.
    """
    # 1. Define constants based on the problem and architecture limitations.
    # We use g=10 for easier representation in the `frac` type.
    # Frac representation: g = 10/1e0 -> {n: 1, d: 1, e: 1}
    g = 10.0
    
    # The lion's speed
    v = 5.0
    
    # The initial distance
    distance = 300.0
    
    # The throwing angle is 60 degrees.
    # We approximate sin(60) and sin(120) to fit `frac` constraints.
    # sin(60) = sqrt(3)/2. Approx sqrt(3) as 26/15.
    # So, sin(60) is approx (26/15)/2 = 13/15. Frac: {n: 13, d: 15, e: 0}
    # sin(120) is the same as sin(60).
    sin_60 = 13 / 15
    sin_120 = 13 / 15

    # 2. Set up the quadratic equation for u: A*u^2 + B*u + C = 0
    # A = sin(2a) = sin(120)
    # B = 2 * v * sin(a) = 2 * v * sin(60)
    # C = -distance * g
    A = sin_120
    B = 2 * v * sin_60
    C = -distance * g

    # 3. Calculate the discriminant: D = B^2 - 4*A*C
    discriminant = B**2 - 4 * A * C

    # 4. Since sqrt() is not available, we implement it with Newton's method.
    def newton_sqrt(n):
        # A good initial guess helps convergence. D is approx 10000.
        x = 100.0
        # 10 iterations provide sufficient precision for this problem.
        for _ in range(10):
            if x == 0: return 0 # Avoid division by zero
            x = 0.5 * (x + n / x)
        return x

    sqrt_of_discriminant = newton_sqrt(discriminant)

    # 5. Solve for u using the quadratic formula (we need the positive speed).
    # u = (-B + sqrt(D)) / (2*A)
    u = (-B + sqrt_of_discriminant) / (2 * A)

    # 6. Calculate memory usage 'z'.
    # The program requires 6 main variables of type 'frac':
    # A, B, C, discriminant, sqrt_of_discriminant, u.
    # Each 'frac' struct (n, d, e) uses 2D + 2D + 2D = 6D of memory.
    # Total memory z = 6 variables * 6D/variable = 36D.
    z = 36

    # 7. Print the components of the final equation as requested.
    print(f"Solving for u in the equation: u = (-B + sqrt(D)) / (2*A)")
    print(f"B = {B:.4f}")
    print(f"D = B^2 - 4*A*C = {discriminant:.4f}")
    print(f"sqrt(D) calculated via Newton's method = {sqrt_of_discriminant:.4f}")
    print(f"A = {A:.4f}")
    print(f"Final calculation: u = ({-B:.4f} + {sqrt_of_discriminant:.4f}) / (2 * {A:.4f})")
    
    # 8. Print the final answer.
    u_rounded = round(u, 1)
    print(f"\nCalculated initial speed u = {u_rounded} m/s")
    print(f"Calculated memory usage z = {z} D")
    print(f"\nFinal Answer (u:z) is: {u_rounded}:{z}")

solve_monkey_and_lion()