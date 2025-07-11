from fractions import Fraction

def solve_wuxing_problem():
    """
    Solves the physics problem based on the Wuxing architecture constraints.
    """

    # Helper function to compute the square root of a Fraction using Newton's method.
    # This simulates a hardware/library function that would be available on Wuxing.
    def sqrt_frac(s_frac, iterations=30):
        # An initial guess for the square root.
        x = s_frac
        half = Fraction(1, 2)
        # Iterate to refine the approximation.
        for _ in range(iterations):
            if x == 0:
                return Fraction(0)
            x = half * (x + s_frac / x)
        return x

    # Step 1: Define physics parameters as Fractions, simulating Wuxing types.
    # int d = 300;     // 5D
    # int v = 5;       // 5D
    # frac g = 98/10e0; // 6D
    d_initial = Fraction(300)
    v_lion = Fraction(5)
    g = Fraction(98, 10)

    # Step 2: Handle trigonometric values. The angle is a=60 degrees.
    # We need sin(60) = sqrt(3)/2. This must be approximated as a rational number.
    # We generate this approximation using our sqrt_frac function.
    # frac sin_a = ...; // 6D
    sqrt_3_approx = sqrt_frac(Fraction(3))
    sin_60 = sqrt_3_approx / Fraction(2)
    
    # In the equation, we also need sin(2a) = sin(120), which equals sin(60).
    sin_120 = sin_60

    # Step 3: Formulate and solve the quadratic equation Au^2 + Bu + C = 0.
    # The equation derived from physics is: u^2*sin(120) + u*(2*v*sin(60)) - d*g = 0
    A = sin_120
    B = Fraction(2) * v_lion * sin_60
    C = -d_initial * g

    # Step 4: Solve for u using the quadratic formula: u = (-B + sqrt(B^2 - 4AC)) / 2A
    discriminant = B**2 - Fraction(4) * A * C
    sqrt_discriminant = sqrt_frac(discriminant)
    
    # We are looking for a positive speed, so we take the '+' root.
    # frac u = ...; // 6D (for the result)
    u_rock = (-B + sqrt_discriminant) / (Fraction(2) * A)

    # Step 5: Calculate memory usage (z).
    # d (int): 5D, v (int): 5D, g (frac): 6D, sin_60 (frac): 6D, u (frac): 6D
    # Total z = 5 + 5 + 6 + 6 + 6 = 28 D
    z = 28
    
    # Step 6: Print the full calculation and final answer.
    u_rounded = round(float(u_rock), 1)
    final_answer_string = f"{u_rounded}:{z}"

    print("Solving u^2*sin(120) + u*(2*v*sin(60)) - d*g = 0")
    print("This is a quadratic equation Au^2 + Bu + C = 0 with:")
    print(f"  d = {float(d_initial)}")
    print(f"  v = {float(v_lion)}")
    print(f"  g = {float(g)}")
    print(f"  sin(60) ≈ {float(sin_60):.6f}")
    print("\nCoefficients:")
    print(f"  A = sin(120) = {float(A):.6f}")
    print(f"  B = 2 * v * sin(60) = {float(B):.6f}")
    print(f"  C = -d * g = {float(C):.6f}")
    
    print("\nUsing the quadratic formula u = (-B + sqrt(B^2 - 4AC)) / 2A:")
    print(f"  Discriminant (B^2 - 4AC) = {float(discriminant):.6f}")
    print(f"  Square root of discriminant ≈ {float(sqrt_discriminant):.6f}")
    print(f"  u = (-({float(B):.6f}) + {float(sqrt_discriminant):.6f}) / (2 * {float(A):.6f})")
    print(f"  u ≈ {float(u_rock):.6f} m/s")

    print("\nFinal Result:")
    print(f"  u (rounded to 0.1) = {u_rounded} m/s")
    print(f"  z (memory usage) = {z} D")
    print(f"\nAnswer in u:z format: {final_answer_string}")

if __name__ == '__main__':
    solve_wuxing_problem()