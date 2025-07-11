def solve_monkey_problem():
    """
    Calculates the initial speed 'u' for the rock to hit the lion,
    adhering to the Wuxing architecture constraints.
    """

    # --- Step 1: Define constants based on the problem ---
    # v: lion's speed (m/s)
    # g: acceleration due to gravity (m/s^2)
    # dist: initial distance to the lion (m)
    v = 5.0
    g = 9.8
    dist = 300.0
    
    # --- Step 2: Handle trigonometric constraint ---
    # The Wuxing architecture has no sin() function. We must use a rational
    # approximation for sin(60) that fits the 'frac' type constraints.
    # sin(60) = sqrt(3)/2. A good approximation is 97/112.
    # frac.n = 97 (fits in signed char)
    # frac.d = 112 (fits in unsigned char)
    sin_60_approx = 97.0 / 112.0

    # --- Step 3: Formulate the quadratic equation Au^2 + Bu + C = 0 ---
    # A = sin(2*60) = sin(120) = sin(60)
    # B = 2 * v * sin(60)
    # C = -dist * g
    A = sin_60_approx
    B = 2 * v * sin_60_approx
    C = -dist * g

    print("The problem is solved using the quadratic equation: A*u^2 + B*u + C = 0")
    print(f"Using rational approximation sin(60) = 97/112 ~= {sin_60_approx:.4f}")
    print(f"Calculated coefficients:")
    print(f"A = sin(60) = {A:.4f}")
    print(f"B = 2 * {v} * sin(60) = {B:.4f}")
    print(f"C = -{dist} * {g} = {C:.1f}")
    print("-" * 20)
    print("Final equation:")
    print(f"{A:.4f} * u^2 + {B:.4f} * u {C:.1f} = 0")
    print("-" * 20)

    # --- Step 4: Handle sqrt constraint ---
    # The Wuxing architecture has no sqrt() function. We implement one using the
    # Babylonian method, which only uses basic arithmetic.
    def my_sqrt(n, iterations=10):
        """Calculates square root using the Babylonian method."""
        # Initial guess can be n/2
        x = n / 2.0
        for _ in range(iterations):
            if x == 0: return 0 # Avoid division by zero
            x = (x + n / x) / 2.0
        return x

    # --- Step 5: Solve the quadratic equation ---
    # u = (-B + sqrt(B^2 - 4AC)) / (2A)
    # We only consider the positive root as speed must be positive.
    discriminant = B**2 - 4 * A * C
    sqrt_discriminant = my_sqrt(discriminant)
    
    print(f"Discriminant (B^2 - 4*A*C) = {discriminant:.4f}")
    print(f"Calculated sqrt(Discriminant) = {sqrt_discriminant:.4f}")
    
    u = (-B + sqrt_discriminant) / (2 * A)
    
    # --- Step 6: Calculate memory usage and format the final answer ---
    # Memory usage 'z' for non-temporary variables (g, v, dist, sin_60_approx)
    # Each 'frac' type variable takes 6D (Decimal positions).
    # z = 4 variables * 6D/variable = 24D
    z = 24
    
    # Final answer format: u:z, with u rounded to one decimal place
    final_answer = f"{u:.1f}:{z}"
    
    print("-" * 20)
    print(f"The required initial speed u is {u:.4f} m/s.")
    print(f"The memory usage z is {z} D.")
    print(f"Final Answer (u:z) = {final_answer}")
    
    # The final output required by the user prompt
    print(f"\n<<<{final_answer}>>>")

solve_monkey_problem()