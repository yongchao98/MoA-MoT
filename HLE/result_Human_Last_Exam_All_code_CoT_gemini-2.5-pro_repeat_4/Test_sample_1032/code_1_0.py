import math

def solve_horse_area():
    """
    Calculates the reachable area for a horse tied to the origin with a taxi-cab rope,
    avoiding a house in the third quadrant.
    """
    
    # Step 1: Define the rope length R
    R = 7.0 / 2.0
    print(f"The taxi-cab rope length is R = 7/2 = {R}")
    
    # Step 2: Calculate the total possible area without obstacles
    total_diamond_area = 2 * R**2
    print(f"The total reachable area without any obstacles is a diamond with area 2 * R^2 = 2 * {R}^2 = {total_diamond_area}")
    
    # Step 3: Calculate the area in the unobstructed quadrants (Q1, Q2, Q4)
    # The area in each quadrant is 1/4 of the total diamond area.
    area_per_quadrant = total_diamond_area / 4.0
    area_q1_q2_q4 = 3 * area_per_quadrant
    print(f"The house is in the third quadrant. Quadrants 1, 2, and 4 are unobstructed.")
    print(f"The area in each of these quadrants is a triangle with area 0.5 * R * R = {area_per_quadrant}")
    print(f"Total area for Q1, Q2, and Q4 = 3 * {area_per_quadrant} = {area_q1_q2_q4}")
    
    # Step 4: Calculate the area in the obstructed third quadrant (Q3)
    print("\nIn Quadrant 3, the house at the origin blocks all direct paths.")
    print("The horse can only reach this area by bending the rope around the corners of the house.")
    
    # The primary pivot points are P1=(-2,0) and P2=(0,-2).
    # Calculate area reachable from pivot P1=(-2,0)
    d_O_P1 = abs(-2 - 0) + abs(0 - 0)
    rem_rope_1 = R - d_O_P1
    print(f"Path 1: Bend at P1=(-2,0). Distance d(O,P1) = {d_O_P1}.")
    print(f"Remaining rope length = {R} - {d_O_P1} = {rem_rope_1}.")
    
    # This forms a new diamond centered at P1 with radius rem_rope_1.
    # We are interested in the semi-diamond in Q3 (where y<0).
    area_from_p1 = (1/2.0) * (2 * rem_rope_1**2)
    print(f"This allows access to a semi-diamond in Q3 with area 0.5 * (2 * {rem_rope_1}^2) = {area_from_p1}.")

    # Calculate area reachable from pivot P2=(0,-2)
    d_O_P2 = abs(0 - 0) + abs(-2 - 0)
    rem_rope_2 = R - d_O_P2
    print(f"Path 2: Bend at P2=(0,-2). Distance d(O,P2) = {d_O_P2}.")
    print(f"Remaining rope length = {R} - {d_O_P2} = {rem_rope_2}.")
    
    # This forms a semi-diamond in Q3 (where x<0).
    area_from_p2 = (1/2.0) * (2 * rem_rope_2**2)
    print(f"This allows access to a semi-diamond in Q3 with area 0.5 * (2 * {rem_rope_2}^2) = {area_from_p2}.")
    
    # Step 5: Check for overlap between the two semi-diamonds in Q3
    print("\nWe need to find the area of the union of these two semi-diamonds.")
    # The region from P1 is defined by |x+2| + |y| <= 1.5, for y<=0.
    # The region from P2 is defined by |x| + |y+2| <= 1.5, for x<=0.
    # A point (x,y) in the intersection must satisfy both inequalities.
    # Let's check a point P(x,y). For it to be in the intersection, we need y >= x+0.5 and y <= x-0.5.
    # Since x+0.5 > x-0.5, there is no value of y that can satisfy both conditions.
    # The intersection is empty.
    area_intersection = 0.0
    print("The two semi-diamonds do not overlap. Their intersection has an area of 0.")
    
    # The total area in Q3 is the sum of the areas of the two non-overlapping regions.
    area_q3 = area_from_p1 + area_from_p2 - area_intersection
    print(f"Total reachable area in Q3 = {area_from_p1} + {area_from_p2} - {area_intersection} = {area_q3}")
    
    # Step 6: Calculate the final total area
    total_area = area_q1_q2_q4 + area_q3
    print("\nCalculating the final answer:")
    print(f"Total Reachable Area = (Area from Q1, Q2, Q4) + (Area from Q3)")
    print(f"Total Reachable Area = {area_q1_q2_q4} + {area_q3} = {total_area}")

solve_horse_area()