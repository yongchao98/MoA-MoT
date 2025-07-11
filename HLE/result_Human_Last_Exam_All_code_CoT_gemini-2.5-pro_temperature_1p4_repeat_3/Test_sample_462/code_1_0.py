import math

def solve_area_between_circles():
    """
    Calculates the area of the region between four circles inscribed in a square.
    """
    # Define the side length of the square
    side_length = 20.0

    # The side length of the square is equal to 4 times the radius of a circle.
    # So, radius = side_length / 4.
    radius = side_length / 4.0

    # Calculate the area of the square
    area_square = side_length ** 2

    # Calculate the total area of the four circles
    area_four_circles = 4 * math.pi * (radius ** 2)

    # Calculate the area of the region between the circles
    area_region = area_square - area_four_circles

    # Print the explanation and the steps
    print(f"Given a square with side length {side_length} cm.")
    print(f"Four circles are inscribed in the corners, touching each other and the sides.")
    print(f"The radius of each circle is the side length divided by 4: {side_length} / 4 = {radius} cm.")
    print("\nTo find the area of the region between the circles, we subtract the area of the four circles from the area of the square.")
    print("\nThe equation is: Area_Region = Area_Square - Area_Four_Circles")
    
    # Print the final equation with the calculated numbers
    print(f"\nFinal Equation with values:")
    print(f"Area_Region = {area_square:.2f} - {area_four_circles:.2f}")
    
    # Print the final result
    final_answer = round(area_region, 2)
    print(f"\nResult: {final_answer} cm^2")

solve_area_between_circles()