import math

def solve_horse_area():
    """
    Calculates the area a horse can reach, tied to a rope at the origin,
    with an L-shaped house as an obstacle.
    """
    
    # The rope has a taxi-cab length of L
    L = 7 / 2
    
    # 1. Calculate the area in the 3 unobstructed quadrants (Q1, Q2, Q4)
    # The area in one quadrant is a triangle with area 0.5 * L * L
    area_one_quadrant = 0.5 * L**2
    area_unobstructed_quads = 3 * area_one_quadrant
    
    # 2. Calculate the area in the obstructed quadrant (Q3)
    
    # The house consists of 3 unit squares, so its area is 3.
    # The house corners are (-2,0),(0,0),(0,-2),(-1,-2),(-1,-1),(-2,-1).
    # It is fully contained within the Q3 taxi-cab circle since the furthest corner
    # (-2,-1) or (-1,-2) has a taxi-cab distance of 3, which is less than L=3.5.
    area_house = 3.0
    
    # The potential area in Q3 without the house blocking the rope
    area_q3_diamond = area_one_quadrant
    
    # a) Directly accessible area in Q3 is the Q3 diamond area minus the house area
    area_q3_direct = area_q3_diamond - area_house
    
    # b) Area accessible by bending the rope around the house's outer corners.
    # The pivot points are P1=(-2, -1) and P2=(-1, -2).
    # The taxi-cab distance from the origin (0,0) to these pivots is the same.
    dist_to_pivot = abs(-2 - 0) + abs(-1 - 0)
    
    # The remaining rope length after reaching a pivot
    L_remaining = L - dist_to_pivot
    
    # From a pivot, the rope can graze in a new taxi-cab circle, but it cannot
    # bend back into the house. This means it can only access 3/4 of the full circle.
    # Area of a full taxi-cab circle is 2 * r^2.
    area_from_one_pivot = (3/4) * (2 * L_remaining**2)
    
    # There are two symmetric pivot points
    area_from_pivots = 2 * area_from_one_pivot
    
    # 3. Calculate the total area
    total_area = area_unobstructed_quads + area_q3_direct + area_from_pivots
    
    # Print the breakdown of the calculation as requested
    print("The total reachable area is calculated by summing disjoint regions:")
    print(f"1. Area in 3 unobstructed quadrants: 3 * (1/2 * {L}^2) = {area_unobstructed_quads}")
    
    print("\nIn Quadrant 3:")
    print(f"2. Directly reachable area = (Area of Q3 diamond) - (Area of house)")
    print(f"   = {area_q3_diamond} - {area_house} = {area_q3_direct}")
    
    print(f"\n3. Area from rope bending around pivots P1(-2,-1) and P2(-1,-2):")
    # We calculate for one pivot and multiply by 2
    dist_p1 = abs(-2) + abs(-1)
    L_rem_p1 = L - dist_p1
    area_p1 = (3/4) * (2 * L_rem_p1**2)
    print(f"   - For pivot P1, remaining rope = {L} - {dist_p1} = {L_rem_p1}")
    print(f"   - Area from P1 = (3/4) * (2 * {L_rem_p1}^2) = {area_p1}")
    print(f"   - Area from P2 is the same: {area_p1}")
    print(f"   Total area from pivots = {area_from_pivots}")
    
    # Print the final equation
    print("\nFinal equation:")
    print(f"Total Area = {area_unobstructed_quads} + {area_q3_direct} + {area_from_pivots}")
    print(f"Total Area = {total_area}")

solve_horse_area()
<<<22.25>>>