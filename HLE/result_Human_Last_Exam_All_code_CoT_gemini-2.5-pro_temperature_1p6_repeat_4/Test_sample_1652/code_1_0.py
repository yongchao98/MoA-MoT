def solve_projectile_problem():
    """
    Calculates the initial speed of a projectile to hit a moving target
    based on the Wuxing architecture constraints.
    """

    # --- Step 1: Define constants and approximations ---
    # According to Wuxing architecture, these would be 'frac' types.
    # We use floating-point numbers here to simulate the logic.
    g = 9.8  # Acceleration due to gravity (m/s^2), as frac 49/5e0
    v = 5.0  # Lion's speed (m/s), as frac 5/1e0
    initial_dist = 300.0  # Initial distance (m), as frac 300/1e0
    angle_deg = 60.0  # Angle in degrees

    # --- Step 2: Implement Wuxing-compatible math functions ---

    def my_sqrt(n, iterations=15):
        """
        Calculates square root using the Babylonian method,
        as the Wuxing 'sqrt' function is not available.
        """
        if n < 0:
            return float('nan') # Not possible in Wuxing, but good practice
        x = n / 2.0
        if x == 0:
            return 0
        for _ in range(iterations):
            x = 0.5 * (x + n / x)
        return x

    # sin(60) = sqrt(3)/2. We must calculate sqrt(3) first.
    # This approximates sqrt(3) using our custom function.
    sqrt_3 = my_sqrt(3.0)
    
    # sin(angle) and sin(2*angle) calculations
    # For a = 60 degrees, sin(a) = sin(60) and sin(2a) = sin(120), and sin(60)=sin(120).
    sin_a = sqrt_3 / 2.0

    # --- Step 3: Formulate and solve the quadratic equation ---
    # The equation for u is: A*u^2 + B*u + C = 0
    # A = sin(2a) / g
    # B = 2 * v * sin(a) / g
    # C = -initial_dist
    
    # In our case, sin(2a) is also sin(120) which is sqrt(3)/2, same as sin(a).
    sin_2a = sin_a

    A = sin_2a / g
    B = (2 * v * sin_a) / g
    C = -initial_dist

    # Calculate the discriminant: D = B^2 - 4*A*C
    discriminant = B**2 - 4 * A * C
    
    # Solve for u using the positive root of the quadratic formula
    # u = (-B + sqrt(D)) / (2*A)
    # We must use our custom sqrt function.
    sqrt_discriminant = my_sqrt(discriminant)
    
    if sqrt_discriminant >= 0:
        u = (-B + sqrt_discriminant) / (2 * A)
    else:
        u = float('nan') # No real solution

    # --- Step 4: Calculate memory usage (z) ---
    # 5 main variables (g, v, initial_dist, sin_a, u)
    # Each 'frac' type costs 6D (2D+2D+2D for n,d,e)
    # z = 5 variables * 6D/variable
    z = 30
    
    # --- Step 5: Print the results ---
    print("Solving for u in the equation: A*u^2 + B*u + C = 0")
    print(f"The equation for the initial speed u is: ({A:.4f})*u^2 + ({B:.4f})*u + ({C:.1f}) = 0")
    print("\nCalculation steps:")
    print(f"B^2 - 4*A*C = ({B:.4f})^2 - 4*({A:.4f})*({C:.1f}) = {discriminant:.4f}")
    print(f"Square root of discriminant = {sqrt_discriminant:.4f}")
    print(f"u = ( -({B:.4f}) + {sqrt_discriminant:.4f} ) / ( 2 * {A:.4f} )")
    print(f"Calculated initial speed u = {u:.4f} m/s")
    
    # Final answer format u:z
    u_rounded = round(u, 1)
    
    print("\n--- Final Answer ---")
    print("u (rounded to 0.1): " + str(u_rounded))
    print("z (memory usage in D): " + str(z))
    print(f"\nFinal result in u:z format: {u_rounded}:{z}")
    
    # The final deliverable for the platform
    final_answer = f"{u_rounded}:{z}"
    return final_answer

final_answer_string = solve_projectile_problem()

# The final output needs to be wrapped in <<<>>>
# print(f"<<<{final_answer_string}>>>")