import math

def solve_for_inner_product():
    """
    Calculates the limiting inner product <b_p, z_p> based on the geometric
    relationship between the three vectors h_p, b_p, and z_p.
    """
    # Given limiting inner products
    c_b = 0.9375  # lim <h_p, b_p>
    c_z = 0.9     # lim <h_p, z_p>

    # The relationship between the cosines of the angles, assuming h_p is
    # a perturbation of b_p towards z_p, is:
    # x = cos(acos(c_b) + acos(c_z))
    # which simplifies to x = c_b*c_z - sqrt(1-c_b^2)*sqrt(1-c_z^2)

    # Calculate the terms of the equation
    term1 = c_b * c_z
    sin_b_sq = 1 - c_b**2
    sin_z_sq = 1 - c_z**2
    
    term2_part1 = math.sqrt(sin_b_sq)
    term2_part2 = math.sqrt(sin_z_sq)
    term2 = term2_part1 * term2_part2
    
    result = term1 - term2

    # Print the equation with the numerical values plugged in, as requested.
    print(f"The desired value is calculated using the formula: c_b * c_z - sqrt(1 - c_b^2) * sqrt(1 - c_z^2)")
    print(f"Value = ({c_b} * {c_z}) - (sqrt(1 - {c_b**2}) * sqrt(1 - {c_z**2}))")
    print(f"Value = {term1} - ({term2_part1} * {term2_part2})")
    print(f"Value = {term1} - {term2}")
    print(f"Final Result: {result}")

solve_for_inner_product()