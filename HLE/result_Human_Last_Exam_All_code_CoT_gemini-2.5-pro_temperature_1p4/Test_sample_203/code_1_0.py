import math

def solve_triangle_median_range():
    """
    Calculates the range of the median m for which angle A is acute,
    based on the given geometric properties.
    """
    # Step 1: Define given constants from the problem
    AD = 12
    AE = 13

    # Step 2: Use Pythagorean theorem in triangle ADE to find DE
    # In our coordinate system A=(0,12), D=(0,0), E=(5,0)
    DE_sq = AE**2 - AD**2
    DE = math.sqrt(DE_sq)
    print(f"The distance DE is sqrt({AE}^2 - {AD}^2) = {DE}")

    # Step 3: Set up the conditions for f, the x-coordinate of the median's foot F
    # From the angle bisector theorem and coordinate geometry, we get relationships
    # that lead to a quadratic equation x^2 - 2fx + (23.8f + 144) = 0 for the
    # coordinates of B and C.

    # Condition 1: For a valid, non-isosceles triangle, the discriminant must be positive.
    # The condition is f^2 - 23.8*f - 144 > 0
    # Let's find the roots of f^2 - 23.8*f - 144 = 0
    # Coefficients for af^2 + bf + c = 0
    a1, b1, c1 = 1, -23.8, -144
    delta1 = b1**2 - 4*a1*c1
    f_root1 = (-b1 - math.sqrt(delta1)) / (2*a1)
    f_root2 = (-b1 + math.sqrt(delta1)) / (2*a1)
    print(f"\nThe condition for a valid triangle leads to f^2 - 23.8f - 144 > 0.")
    print(f"The roots of this quadratic are f = {f_root1} and f = {f_root2}.")
    print(f"So, f must be in the range (-inf, {f_root1}) U ({f_root2}, inf).")

    # Condition 2: For angle A to be acute, we have bc + 144 > 0.
    # This leads to 23.8f + 288 > 0
    f_lower_bound_for_acute_A = -288 / 23.8
    f_lower_bound_num, f_lower_bound_den = -1440, 119
    print(f"\nThe condition for angle A to be acute leads to f > -288 / 23.8.")
    print(f"This simplifies to f > {f_lower_bound_num}/{f_lower_bound_den} (approx {f_lower_bound_for_acute_A:.4f}).")

    # Step 4: Combine the conditions for f and find the corresponding range for m.
    # m = sqrt(f^2 + AD^2)

    # First valid interval for f: f > f_root2 (f > 28.8)
    # This interval satisfies f > -1440/119.
    m_lower1 = math.sqrt(f_root2**2 + AD**2)
    print(f"\nFor the first range f > {f_root2}, the corresponding range for m is m > sqrt({f_root2}^2 + {AD}^2).")
    print(f"m > {m_lower1}")

    # Second valid interval for f: f_lower_bound_for_acute_A < f < f_root1
    # -1440/119 < f < -5
    # The function m(f) is decreasing for f < 0.
    # So, m(f_root1) < m < m(f_lower_bound_for_acute_A)
    m_lower2 = math.sqrt(f_root1**2 + AD**2)
    m_upper2_sq = (f_lower_bound_num / f_lower_bound_den)**2 + AD**2
    m_upper2 = math.sqrt(m_upper2_sq)
    
    m_upper2_num = 12 * 169
    m_upper2_den = 119

    print(f"\nFor the second range {f_lower_bound_num}/{f_lower_bound_den} < f < {f_root1}, the corresponding range for m is:")
    print(f"sqrt({f_root1}^2 + {AD}^2) < m < sqrt(({f_lower_bound_num}/{f_lower_bound_den})^2 + {AD}^2)")
    print(f"{m_lower2} < m < {m_upper2_num}/{m_upper2_den} (approx {m_upper2:.4f})")

    # Final result summary
    print("\n--- Summary ---")
    print(f"The range of values for m for which angle A will be acute is the union of two intervals:")
    print(f"({m_lower2}, {m_upper2_num}/{m_upper2_den}) U ({m_lower1}, infinity)")

solve_triangle_median_range()