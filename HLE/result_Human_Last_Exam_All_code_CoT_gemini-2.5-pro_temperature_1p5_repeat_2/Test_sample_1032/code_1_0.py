def solve_horse_area():
    """
    Calculates the reachable area for the horse.
    """
    # Length of the rope
    L = 3.5

    # --- Step 1: Area in Quadrants 1, 2, and 4 ---
    # Area in one quadrant is a triangle with base L and height L.
    # Area = 0.5 * L * L. For three quadrants:
    area_q124 = 3 * (0.5 * L**2)
    print(f"Step 1: The unobstructed area in Quadrants 1, 2, and 4 is 3 * (1/2 * {L}^2) = {area_q124}")

    # --- Step 2: Directly reachable area in Quadrant 3 ---
    # This is the area of the diamond |x|+|y|<=L within the region (x>=-1 or y>=-1) in Q3,
    # minus the area of the house.
    # Area of diamond in strip [-1,0]x(inf,0) = Integral from -1 to 0 of (x+3.5) dx = 3
    # Area of diamond in strip (inf,0)x[-1,0] = Integral from -1 to 0 of (y+3.5) dy = 3
    # Area of diamond in square [-1,0]x[-1,0] = 1
    # Area of union = 3 + 3 - 1 = 5
    area_q3_unblocked_gross = 5.0
    # The house has an area of 3 unit squares.
    house_area = 3.0
    area_q3_unblocked_net = area_q3_unblocked_gross - house_area
    print(f"Step 2: The directly reachable area in Quadrant 3 is (5 - 3) = {area_q3_unblocked_net}")

    # --- Step 3 & 4: Area reachable by bending the rope (and overlap) ---
    # This is the area in the "shadow" region (x<-1 and y<-1)
    
    # Area reachable from pivot A(-2,0)
    # This area is a triangle with vertices (-2,-1), (-2.5,-1), (-2,-1.5)
    area_A_reach = 0.5 * 0.5 * 0.5
    print(f"Step 3a: The extra area reachable by pivoting around A(-2,0) is {area_A_reach}")

    # Area reachable from pivot C(0,-2)
    # This area is composed of two triangles.
    # T1 vertices: (-1, -1.5), (-1.5, -2), (-1, -2) -> Area = 0.125
    # T2 vertices: (-1.5, -2), (-1, -2.5), (-1, -2) -> Area = 0.125
    area_C_reach = 0.125 + 0.125
    print(f"Step 3b: The extra area reachable by pivoting around C(0,-2) is {area_C_reach}")
    
    # Overlap between the two regions above
    # This area is a triangle with vertices (-2,-1.5), (-1.5,-2), (-2,-2)
    area_overlap = 0.5 * 0.5 * 0.5
    print(f"Step 4: The overlap between the two pivoted regions is {area_overlap}")

    # --- Step 5: Final Summation ---
    total_area = area_q124 + area_q3_unblocked_net + area_A_reach + area_C_reach - area_overlap
    print("\n--- Final Calculation ---")
    print(f"Total Area = (Area Q1,2,4) + (Direct Q3 Area) + (Area via A) + (Area via C) - (Overlap)")
    print(f"Total Area = {area_q124} + {area_q3_unblocked_net} + {area_A_reach} + {area_C_reach} - {area_overlap}")
    print(f"Total Area = {total_area}")
    return total_area

solve_horse_area()
<<<20.625>>>