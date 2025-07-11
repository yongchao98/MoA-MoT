import math

# Plan:
# 1. Set up the physics problem which results in a quadratic equation for the unknown initial speed 'u'.
#    The equation is: A*u^2 + B*u + C = 0, where:
#    A = sin(2*a)
#    B = 2*v*sin(a)
#    C = -dist*g
# 2. Since the Wuxing architecture lacks sqrt() and sin() functions, we will simulate the solution.
#    - Use a pre-calculated fractional value for sin(60) and sin(120).
#    - Implement a custom sqrt() function using Newton's method.
# 3. Solve for 'u' using the quadratic formula with our custom sqrt function.
# 4. Calculate the memory usage 'z' based on the specified Wuxing data types.
# 5. Print the equation with all numerical values, then output the final answer in the required format.

def solve_projectile_velocity():
    """
    Calculates the initial velocity 'u' for the rock to hit the lion.
    """

    # Custom square root function to comply with architectural constraints (no math.sqrt).
    def custom_sqrt(number, iterations=15):
        """Calculates square root using Newton's method."""
        if number < 0:
            # Should not happen for a real solution in this physics problem
            return 0 
        
        # Start with an initial guess
        x = number / 2.0 if number > 0 else 0.0
        
        # Iterate to improve precision
        for _ in range(iterations):
            if x == 0.0:
                break
            x = (x + number / x) / 2.0
        return x

    # --- Step 1: Define problem constants ---
    # These variables would be 'int' or 'frac' types on the Wuxing system.
    v = 5.0            # Lion's speed in m/s (can be stored as Wuxing int)
    dist = 300.0       # Initial distance in m (can be stored as Wuxing int)
    g = 9.8            # Acceleration due to gravity (must be Wuxing frac, e.g., 98/10e0)
    angle_a = 60       # Angle of throw in degrees
    
    # Pre-calculated sine values, as sin() is not available.
    # sin(60) = sin(120) is approx 0.866. This is stored as a Wuxing 'frac'.
    sin_a = 0.866
    sin_2a = 0.866

    # --- Step 2: Calculate coefficients for the quadratic equation A*u^2 + B*u + C = 0 ---
    A = sin_2a
    B = 2 * v * sin_a
    C = -dist * g

    # --- Step 3: Solve the quadratic equation for u ---
    # We need the positive solution for a forward throw.
    # u = (-B + sqrt(B^2 - 4*A*C)) / (2*A)
    discriminant = B**2 - 4 * A * C
    
    if discriminant < 0:
        print("The rock can never hit the lion under these conditions.")
        return

    sqrt_of_discriminant = custom_sqrt(discriminant)
    
    u = (-B + sqrt_of_discriminant) / (2 * A)

    # --- Step 4: Output the equation with numerical values ---
    # As per instructions, show the calculation with the numbers plugged in.
    print("Solving the quadratic equation: u = (-B + sqrt(B^2 - 4*A*C)) / (2*A)")
    print(f"u = (-{B:.3f} + sqrt({B:.3f}^2 - 4 * {A:.3f} * {C:.1f})) / (2 * {A:.3f})")

    # --- Step 5: Calculate memory usage (z) ---
    # Data types in Wuxing: int (5D), char (2D), frac (6D)
    # v (as int): 5D
    # dist (as int): 5D
    # g (as frac): 6D
    # sin_a (as frac): 6D
    # u (result, as frac): 6D
    z = 5 + 5 + 6 + 6 + 6

    # --- Step 6: Format final answer ---
    u_rounded = round(u, 1)
    
    print(f"<<<{u_rounded}:{z}>>>")

solve_projectile_velocity()