import math

def solve_horse_area():
    """
    Calculates the reachable area for a horse tied with a rope of taxi-cab length.
    """
    
    # Rope length
    L = 7.0 / 2.0

    # Step 1: Calculate the unobstructed area in Quadrants 1, 2, and 4.
    # The area in each quadrant is a triangle with base L and height L.
    # Area of one such triangle = 1/2 * L * L
    area_one_quadrant = 0.5 * L**2
    area_q1_q2_q4 = 3 * area_one_quadrant

    print(f"The rope length is L = {L}.")
    print("\nStep 1: Calculate the area in Quadrants 1, 2, and 4 (unobstructed).")
    print(f"The area in each of these three quadrants is a triangle with area 1/2 * L^2.")
    print(f"Area = 3 * (1/2 * {L}^2) = 3 * {area_one_quadrant} = {area_q1_q2_q4}")

    # Step 2: Calculate the area in Quadrant 3, considering the house obstacle.
    # The rope pivots around the corners of the house.

    print("\nStep 2: Calculate the area in Quadrant 3 (obstructed).")
    print("The rope pivots around the corners of the house. We calculate the area accessible from each pivot and handle overlaps.")

    # Distance to primary pivots A(-2,0), C(0,-2), and E(-1,-1)
    dist_to_A = abs(-2) + abs(0)  # Pivot A(-2,0)
    dist_to_C = abs(0) + abs(-2)  # Pivot C(0,-2)
    # Path to E(-1,-1) goes along the boundary: |-1|+|0| + |0|+|-1| = 2
    dist_to_E = abs(-1) + abs(-1) 

    # Remaining rope length from these pivots
    rem_L_AC = L - dist_to_A
    rem_L_E = L - dist_to_E
    
    # Area from pivot A(-2,0)
    area_from_A = 0.5 * rem_L_AC**2
    print(f"\nArea from pivot A(-2,0):")
    print(f"  - Distance from origin to A is {dist_to_A}.")
    print(f"  - Remaining rope length is {L} - {dist_to_A} = {rem_L_AC}.")
    print(f"  - New area is a triangle: 1/2 * {rem_L_AC}^2 = {area_from_A}")

    # Area from pivot C(0,-2)
    area_from_C = 0.5 * rem_L_AC**2
    print(f"\nArea from pivot C(0,-2):")
    print(f"  - Distance from origin to C is {dist_to_C}.")
    print(f"  - Remaining rope length is {L} - {dist_to_C} = {rem_L_AC}.")
    print(f"  - New area is a triangle: 1/2 * {rem_L_AC}^2 = {area_from_C}")

    # Area from pivot E(-1,-1)
    area_from_E = 0.5 * rem_L_E**2
    print(f"\nArea from pivot E(-1,-1):")
    print(f"  - Distance from origin to E (around the house) is {dist_to_E}.")
    print(f"  - Remaining rope length is {L} - {dist_to_E} = {rem_L_E}.")
    print(f"  - New area is a triangle: 1/2 * {rem_L_E}^2 = {area_from_E}")

    # Step 3: Handle overlaps using inclusion-exclusion.
    # The overlap between the regions from A and E is the region accessible from pivot F(-2,-1).
    # The overlap between the regions from C and E is the region accessible from pivot D(-1,-2).
    # The overlap between regions from A and C is empty.
    
    dist_to_F = dist_to_A + 1  # Path O -> A -> F(-2,-1)
    rem_L_F = L - dist_to_F
    area_overlap_AE = 0.5 * rem_L_F**2

    dist_to_D = dist_to_C + 1  # Path O -> C -> D(-1,-2)
    rem_L_D = L - dist_to_D
    area_overlap_CE = 0.5 * rem_L_D**2
    
    print("\nStep 3: Handle overlaps.")
    print("The overlap between areas from A and E corresponds to the area from pivot F(-2,-1).")
    print(f"  - Distance to F is {dist_to_F}, remaining rope is {rem_L_F}.")
    print(f"  - Overlap area = 1/2 * {rem_L_F}^2 = {area_overlap_AE}")

    print("The overlap between areas from C and E corresponds to the area from pivot D(-1,-2).")
    print(f"  - Distance to D is {dist_to_D}, remaining rope is {rem_L_D}.")
    print(f"  - Overlap area = 1/2 * {rem_L_D}^2 = {area_overlap_CE}")
    
    # Area in Q3 is Area(A) + Area(C) + Area(E) - Area(A and E) - Area(C and E)
    area_q3 = area_from_A + area_from_C + area_from_E - area_overlap_AE - area_overlap_CE
    
    print("\nTotal area in Quadrant 3 is calculated using the Principle of Inclusion-Exclusion:")
    print(f"Area(Q3) = Area(A) + Area(C) + Area(E) - Overlap(A,E) - Overlap(C,E)")
    print(f"Area(Q3) = {area_from_A} + {area_from_C} + {area_from_E} - {area_overlap_AE} - {area_overlap_CE} = {area_q3}")

    # Step 4: Final Calculation
    total_area = area_q1_q2_q4 + area_q3
    
    print("\nStep 4: Final Calculation.")
    print("Total Area = Area(Q1,Q2,Q4) + Area(Q3)")
    print(f"Total Area = {area_q1_q2_q4} + {area_q3} = {total_area}")

solve_horse_area()
<<<21.5>>>