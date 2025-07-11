import math

def solve_cube_area_ratio():
    """
    Calculates and prints the exact ratio of the specified area D to the
    total surface area S of a cube.

    The final ratio is (2π + 3√3 - 3) / 12.
    """

    # Define the components of the final formula for the numerator
    coefficient_pi = 2
    coefficient_sqrt3 = 3
    constant_term = -3

    # Define the denominator
    denominator = 12

    # Since the user wants to see the numbers in the equation,
    # we format it as a string for display.
    numerator_str = f"({coefficient_pi}*π + {coefficient_sqrt3}*√3 + ({constant_term}))"

    print("The area of D divided by the area of S is given by the exact expression:")
    print(f"Ratio = {numerator_str} / {denominator}")

    # For context, we can also print the numerical approximation.
    numerator_value = coefficient_pi * math.pi + coefficient_sqrt3 * math.sqrt(3) + constant_term
    ratio_value = numerator_value / denominator
    print("\nThe approximate numerical value of this ratio is:")
    print(ratio_value)

solve_cube_area_ratio()