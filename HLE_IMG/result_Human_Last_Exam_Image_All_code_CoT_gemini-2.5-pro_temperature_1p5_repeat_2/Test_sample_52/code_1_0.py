import math

def calculate_and_print_periods():
    """
    Calculates and prints the tiling periods for the four specified shapes.

    The period is the area of the smallest repeating unit cell, measured in
    units of the area of a single hexagon.
    """

    print("Calculating the tiling periods for four shapes on a hexagonal grid.")
    print("The area of a single hexagon is considered 1 unit.\n")

    # --- Shape 1: 13, 31, 23 ---
    # This shape is a triangle connecting the centers of three adjacent hexagons.
    # Its area is 1/2 of a hexagon's area.
    # A unit cell that can tile the plane is a rhombus made of 2 such tiles.
    area1 = 0.5
    num_tiles1 = 2
    period1 = num_tiles1 * area1
    print("--- Shape 1: 13, 31, 23 ---")
    print("This tile is a triangle with an area of 1/2 hexagon units.")
    print("The smallest repeating unit cell is a rhombus made of 2 tiles.")
    print(f"The defining equation for the unit cell area is: {num_tiles1} * {area1} = {int(period1)}")
    print(f"The period is {int(period1)}.\n")


    # --- Shape 2: 10, 4, 23, 31 ---
    # This shape is interpreted as a rectangle with area equal to 1 hexagon unit.
    # A rectangle can tile the plane by translation, so it is its own unit cell.
    area2 = 1.0
    num_tiles2 = 1
    period2 = num_tiles2 * area2
    print("--- Shape 2: 10, 4, 23, 31 ---")
    print("This tile is a rectangle with an area of 1 hexagon unit.")
    print("It can tile the plane by translation alone, so the unit cell is the tile itself.")
    print(f"The defining equation for the unit cell area is: {num_tiles2} * {area2} = {int(period2)}")
    print(f"The period is {int(period2)}.\n")


    # --- Shape 3: 5, 15, 17, 19, 21, 7 ---
    # This shape is a hexagon with an area of 3/4 of a hexagon unit.
    # We need to find the smallest integers n, m such that n * (3/4) = m.
    # This is satisfied for n=4, which gives m=3.
    area3 = 0.75
    num_tiles3 = 4
    period3 = num_tiles3 * area3
    print("--- Shape 3: 5, 15, 17, 19, 21, 7 ---")
    print("This tile is a hexagon with an area of 3/4 hexagon units.")
    print("The smallest repeating unit cell is made of 4 tiles.")
    print(f"The defining equation for the unit cell area is: {num_tiles3} * {area3} = {int(period3)}")
    print(f"The period is {int(period3)}.\n")


    # --- Shape 4: 4, 5, 14, 23, 18, 19, 21, 22, 31, 30, 9, 8, 13 ---
    # This shape is the union of the three hexagons, with an area of 3 hexagon units.
    # The smallest unit cell that tiles by translation is made of 2 such tiles (with rotation).
    area4 = 3.0
    num_tiles4 = 2
    period4 = num_tiles4 * area4
    print("--- Shape 4: a 13-point polygon ---")
    print("This tile is a Y-shape with an area of 3 hexagon units.")
    print("The smallest repeating unit cell is a large hexagon made of 2 tiles (with rotation).")
    print(f"The defining equation for the unit cell area is: {num_tiles4} * {area4} = {int(period4)}")
    print(f"The period is {int(period4)}.\n")


    # --- Final Answer ---
    final_answer = f"{int(period1)},{int(period2)},{int(period3)},{int(period4)}"
    print("The four periods separated by commas are:")
    print(final_answer)

# Execute the function
calculate_and_print_periods()