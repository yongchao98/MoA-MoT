import math

def solve_horse_area():
    """
    Calculates the total area a horse can reach on a rope with a given
    taxi-cab length, constrained by a house-shaped obstacle.
    """
    # The length of the rope
    rope_length = 7.0 / 2.0

    # The area of the house is 3 unit squares.
    house_area = 3.0

    # --- Calculation ---

    # 1. Area of the initial diamond region, divided into quadrants.
    # The area of a full diamond is 2 * L^2. A single quadrant is 0.5 * L^2.
    area_per_quadrant = 0.5 * rope_length**2
    
    # Quadrants 1, 2, and 4 are unobstructed.
    area_unobstructed_quads = 3 * area_per_quadrant

    # 2. In Quadrant 3, subtract the house area.
    area_q3_initial = area_per_quadrant - house_area
    
    # 3. Calculate extra areas from rope wrapping around the house's outer corners.
    # The pivot corners are (-2,0), (0,-2), (-2,-1), and (-1,-2).

    # Extra area from pivot at P1=(-2, 0).
    # Taxi-cab distance from origin to P1 is 2.
    rope_rem_1 = rope_length - 2.0
    extra_area_1 = 0.5 * rope_rem_1**2

    # Extra area from pivot at P2=(0, -2).
    # Taxi-cab distance from origin to P2 is 2.
    rope_rem_2 = rope_length - 2.0
    extra_area_2 = 0.5 * rope_rem_2**2
    
    # Extra area from pivot at P3=(-2, -1).
    # Shortest path is origin -> (-2,0) -> (-2,-1). Distance = 2 + 1 = 3.
    rope_rem_3 = rope_length - 3.0
    extra_area_3 = 0.5 * rope_rem_3**2
    
    # Extra area from pivot at P4=(-1, -2).
    # Shortest path is origin -> (0,-2) -> (-1,-2). Distance = 2 + 1 = 3.
    rope_rem_4 = rope_length - 3.0
    extra_area_4 = 0.5 * rope_rem_4**2

    # 4. Sum all the calculated areas.
    total_area = (
        area_unobstructed_quads + 
        area_q3_initial + 
        extra_area_1 + 
        extra_area_2 + 
        extra_area_3 + 
        extra_area_4
    )
    
    # --- Output the final result as a single equation ---
    
    print("The total reachable area is the sum of the areas in the unobstructed quadrants, the obstructed quadrant, and the extra areas from the rope wrapping around the corners.")
    print("Final Equation:")
    print(f"Area = (Area Q1,2,4) + (Area Q3 initial) + (Extra Area 1) + (Extra Area 2) + (Extra Area 3) + (Extra Area 4)")
    # Using 'g' format to remove trailing zeros for cleaner display
    print(f"Area = {area_unobstructed_quads:g} + {area_q3_initial:g} + {extra_area_1:g} + {extra_area_2:g} + {extra_area_3:g} + {extra_area_4:g} = {total_area:g}")

solve_horse_area()