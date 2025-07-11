import math

def solve_horse_area():
    """
    Calculates the area a horse can reach, tied to a rope with a specific
    taxi-cab length, avoiding a house-shaped obstacle.
    """
    # Step 1: Define constants
    rope_length = 3.5  # L = 7/2
    house_area = 3.0

    # Step 2: Calculate the area in the three unobstructed quadrants (Q1, Q2, Q4)
    # The area in one quadrant is a triangle with base=L and height=L.
    # Area = 1/2 * base * height = 1/2 * L^2
    unobstructed_quadrant_area = 0.5 * rope_length**2
    area_q1_q2_q4 = 3 * unobstructed_quadrant_area

    # Step 3: Calculate the area in the third quadrant (Q3)
    # 3a. Directly reachable area in Q3
    # This is the area of the Q3 triangle minus the house area.
    # The entire house is within the |x|+|y|<=3.5 diamond.
    direct_area_q3 = unobstructed_quadrant_area - house_area

    # 3b. Area reachable by bending the rope around the house corners
    # The corners act as new pivots.
    # Pivot P1=(-2,0), distance from origin is 2.
    # Pivot P2=(0,-2), distance from origin is 2.
    # Pivot P3=(-1,-1), path length from origin is 2.
    dist_to_pivots = 2.0
    remaining_rope = rope_length - dist_to_pivots # 3.5 - 2.0 = 1.5

    # The additional area from each pivot is a right triangle with side length `remaining_rope`.
    # Its area is 1/2 * base * height = 1/2 * remaining_rope^2
    area_from_one_pivot = 0.5 * remaining_rope**2

    # We must account for overlaps between the areas from the pivots.
    # Overlap between P1 and P3 areas is a small triangle.
    # Its vertices are (-2,-1), (-2.5,-1), (-2,-1.5). Side lengths are 0.5.
    # Overlap between P2 and P3 areas is a similar triangle with side lengths 0.5.
    # Overlap between P1 and P2 areas is zero.
    overlap_side_length = 0.5
    area_of_one_overlap = 0.5 * overlap_side_length**2

    # Total bend area = Area(B1)+Area(B2)+Area(B3) - Area(B1nB3) - Area(B2nB3)
    total_bend_area = 3 * area_from_one_pivot - 2 * area_of_one_overlap

    # Total area in Q3 is the sum of direct and bend areas
    total_area_q3 = direct_area_q3 + total_bend_area

    # Step 4: Calculate total area
    total_area = area_q1_q2_q4 + total_area_q3

    # Step 5: Print the final equation with all numbers
    # The calculation is: (Area Q1,2,4) + (Direct Area Q3) + (Bend Area Q3)
    # Area Q1,2,4 = 3 * (1/2 * 3.5^2)
    # Direct Area Q3 = (1/2 * 3.5^2) - 3
    # Bend Area Q3 = 3 * (1/2 * 1.5^2) - 2 * (1/2 * 0.5^2)
    
    print("The final calculation is composed of three parts:")
    print("1. Area of three unobstructed quadrants (Q1, Q2, Q4)")
    print("2. Directly reachable area in Q3 (unobstructed part of Q3 minus the house)")
    print("3. Extra area in Q3 reachable by bending the rope around the house corners")
    print("\nFinal Equation:")
    print(f"Total Area = ({area_q1_q2_q4}) + ({direct_area_q3}) + ({total_bend_area})")
    print(f"Total Area = {area_q1_q2_q4 + direct_area_q3 + total_bend_area}")

solve_horse_area()