import math

def solve_horse_area():
    """
    Calculates the area a horse can reach, tied to the origin by a rope of a given
    taxi-cab length, with a house-shaped obstacle.
    """

    # 1. Define the parameters
    rope_length = 7.0 / 2.0  # L = 3.5

    # 2. Calculate the total reachable area without obstacles
    # The area of a taxi-cab circle |x| + |y| <= L is 2 * L^2.
    total_area = 2 * rope_length**2

    # 3. Calculate the area of the house (obstacle)
    # The house is described by vertices:
    # (-2, 0), (0, 0), (0, -2), (-1, -2), (-1, -1), (-2, -1)
    # This shape can be decomposed into two rectangles in the third quadrant:
    # Rectangle 1: from x=-2 to x=0, and y=-1 to y=0.
    # Width = 0 - (-2) = 2. Height = 0 - (-1) = 1. Area = 2 * 1 = 2.
    # Rectangle 2: from x=-1 to x=0, and y=-2 to y=-1.
    # Width = 0 - (-1) = 1. Height = -1 - (-2) = 1. Area = 1 * 1 = 1.
    house_area_rect1 = 2.0
    house_area_rect2 = 1.0
    house_area = house_area_rect1 + house_area_rect2

    # 4. Explain the key insight
    # The shortest path from the origin (a corner of the house) around the
    # axis-aligned house to any point (x, y) in the third quadrant has a
    # taxi-cab length of |x| + |y|. This means the house doesn't create any
    # "shadow" region beyond its own footprint. The reachable region is simply
    # the full taxi-cab circle minus the area of the house itself.
    
    # 5. Calculate the final reachable area
    reachable_area = total_area - house_area
    
    # 6. Print the step-by-step calculation
    print("Step 1: Calculate the total possible area.")
    print(f"The rope has a taxi-cab length of L = {rope_length}.")
    print(f"The area of the taxi-cab circle |x| + |y| <= {rope_length} is 2 * L^2.")
    print(f"Total Area = 2 * {rope_length}^2 = {total_area}")
    print("\nStep 2: Calculate the area of the house (the obstacle).")
    print(f"The house consists of two rectangles with areas {house_area_rect1} and {house_area_rect2}.")
    print(f"House Area = {house_area_rect1} + {house_area_rect2} = {house_area}")
    print("\nStep 3: Calculate the final reachable area.")
    print("The final area is the total area minus the house area.")
    print(f"Reachable Area = {total_area} - {house_area} = {reachable_area}")

solve_horse_area()
<<<21.5>>>