import math

def solve_horse_area():
    """
    Calculates the area a horse can reach, tied with a rope of taxi-cab length,
    with an L-shaped house as an obstacle.
    """
    # The rope has a taxi-cab length of 7/2.
    rope_length = 7.0 / 2.0

    # Step 1: Calculate the unobstructed area.
    # The region is |x| + |y| <= L, which is a diamond shape.
    # Its area is 2 * L^2.
    area_diamond = 2 * rope_length**2

    # Step 2: Define the house area to be subtracted.
    # The house consists of three 1x1 unit squares. Its total area is 3.
    # This obstacle lies entirely within the main diamond.
    area_house = 3.0

    # Step 3: Calculate extra area from the rope bending around corners.
    # The corners act as new pivot points for the rope.

    # Pivot around corner C1 = (-2, 0).
    # Taxi-cab distance from origin to C1 is 2.0.
    remaining_rope_c1 = rope_length - 2.0
    # The extra area is a quarter-diamond, with area 0.5 * r^2.
    area_add_c1 = 0.5 * remaining_rope_c1**2

    # Pivot around corner C2 = (0, -2). Symmetric to C1.
    # Taxi-cab distance from origin to C2 is 2.0.
    remaining_rope_c2 = rope_length - 2.0
    area_add_c2 = 0.5 * remaining_rope_c2**2

    # Pivot around corner C3 = (-2, -1), reached via C1.
    # Path length to C3 is d(O, C1) + d(C1, C3) = 2.0 + 1.0 = 3.0.
    remaining_rope_c3 = rope_length - 3.0
    area_add_c3 = 0.5 * remaining_rope_c3**2

    # Pivot around corner C4 = (-1, -2), reached via C2.
    # Path length to C4 is d(O, C2) + d(C2, C4) = 2.0 + 1.0 = 3.0.
    remaining_rope_c4 = rope_length - 3.0
    area_add_c4 = 0.5 * remaining_rope_c4**2
    
    # Step 4: Sum all the parts to get the total area.
    # The added areas are disjoint, so we can simply sum them.
    total_area = area_diamond - area_house + area_add_c1 + area_add_c2 + area_add_c3 + area_add_c4

    # Print the final equation with the calculated numbers
    print("The total area is calculated by the following equation:")
    print(f"Total Area = (Initial Diamond Area) - (House Area) + (Area from Corner 1) + (Area from Corner 2) + (Area from Corner 3) + (Area from Corner 4)")
    print(f"Total Area = {area_diamond} - {area_house} + {area_add_c1} + {area_add_c2} + {area_add_c3} + {area_add_c4}")
    print("\nThe final calculated area is:")
    print(total_area)

solve_horse_area()
<<<24.0>>>