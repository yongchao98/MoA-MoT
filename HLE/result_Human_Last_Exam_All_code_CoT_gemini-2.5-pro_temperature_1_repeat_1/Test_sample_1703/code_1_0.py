import math

def calculate_triangle_area_function():
    """
    This function calculates the area of the triangle T(t) as a function of time t.
    """
    # Given parameters
    R = 10  # Radius of the circle
    v = 1   # Speed of the vertices

    # --- Step 1: Explain the setup ---
    print("Step 1: Simplify the problem")
    print("The area of the triangle is independent of the hexagon's rotation.")
    print("We can analyze the geometry in a stationary frame where the hexagon does not rotate.")
    print("-" * 30)

    # --- Step 2: Calculate the initial area (the constant term) ---
    print("Step 2: Calculate the initial area at t=0 (A_0)")
    # Side length of the hexagon is equal to the radius of the circumscribed circle.
    s = R
    print(f"The hexagon's side length, s = R = {s}")

    # The vertices of T(0) are midpoints of alternating sides of H.
    # This forms an equilateral triangle with side length a = (3/2)s.
    a = (3/2) * s
    print(f"The initial triangle T(0) is equilateral with side length a = (3/2)s = {a}")

    # The area of an equilateral triangle is (sqrt(3)/4) * a^2.
    # This is the constant term in our area function.
    # A_0 = (sqrt(3)/4) * ((3/2)s)^2 = (sqrt(3)/4) * (9/4)s^2 = 9*sqrt(3)/16 * s^2
    A0_coeff_numerator = 9 * (s**2)
    A0_coeff_denominator = 16
    A0_numerical = (9 * math.sqrt(3) / 16) * (s**2)

    print(f"The initial area A_0 = (9 * sqrt(3) / 16) * s^2")
    print(f"A_0 = ({A0_coeff_numerator} * sqrt(3)) / {A0_coeff_denominator}")
    print(f"Numerically, A_0 ≈ {A0_numerical:.4f}")
    print("-" * 30)

    # --- Step 3: Calculate the time-dependent term ---
    print("Step 3: Calculate the coefficient of the t^2 term (C)")
    # Through coordinate geometry (as detailed in the thinking process),
    # the area function is found to be A(t) = A_0 + C * t^2.
    # The coefficient C is (3*sqrt(3)/4) * v^2.
    C_coeff_numerator = 3 * (v**2)
    C_coeff_denominator = 4
    C_numerical = (3 * math.sqrt(3) / 4) * (v**2)
    
    print(f"The area increases over time. The change is proportional to t^2.")
    print(f"The coefficient C of the t^2 term is (3 * sqrt(3) / 4) * v^2")
    print(f"C = ({C_coeff_numerator} * sqrt(3)) / {C_coeff_denominator}")
    print(f"Numerically, C ≈ {C_numerical:.4f}")
    print("-" * 30)

    # --- Step 4: Present the final formula ---
    print("Step 4: Final Area Function A(t)")
    print("Combining the constant and time-dependent parts, the area of the triangle T(t) is:")
    
    # Using integer arithmetic for exact fractions
    final_A0_num = 9 * R**2
    final_A0_den = 16
    common_divisor_A0 = math.gcd(final_A0_num, final_A0_den)
    final_A0_num //= common_divisor_A0
    final_A0_den //= common_divisor_A0

    final_C_num = 3 * v**2
    final_C_den = 4

    print(f"\nArea(t) = ({final_A0_num}√3 / {final_A0_den}) + ({final_C_num}√3 / {final_C_den}) * t²")
    print(f"Numerically, this is approximately: Area(t) ≈ {A0_numerical:.4f} + {C_numerical:.4f} * t²")

calculate_triangle_area_function()