import math

def solve_horse_area():
    """
    Calculates the reachable area for the horse with a rope of taxi-cab length,
    considering a house as an obstacle.
    """
    # Given rope taxi-cab length
    rope_length = 7.0 / 2.0

    # Step 1: Calculate the total area reachable without any obstacles.
    # The area of the region |x| + |y| <= L is 2 * L^2.
    initial_area = 2 * rope_length**2

    # Step 2: Define the area of the house.
    # The house is an L-shaped region with an area equivalent to 3 unit squares.
    house_area = 3.0

    # Step 3: Calculate the area of the "shadow" regions.
    # The relevant pivot corners are C1=(-2,0), C2=(0,-2), and C3=(-1,-1).
    # The taxi-cab distance from the origin (0,0) to each of these corners is 2.
    dist_to_pivot = 2.0
    remaining_rope = rope_length - dist_to_pivot

    # Each pivot creates a quarter-diamond shaped shadow region with area 0.5 * r^2.
    # Since the remaining rope length is the same for all three pivots, their shadow areas are equal.
    shadow_area_each = 0.5 * remaining_rope**2
    
    # Sum of the three individual shadow areas (A1, A2, A3)
    sum_of_shadows = 3 * shadow_area_each

    # Step 4: Use the Principle of Inclusion-Exclusion for the total shadow area.
    # The overlap between any two regions is a triangle with legs of length 0.5.
    # Overlap area = 0.5 * (0.5)^2 = 0.125.
    overlap_area = 0.125
    
    # There are 3 pairwise overlaps.
    sum_of_pairwise_overlaps = 3 * overlap_area
    
    # The triple overlap region is the same as any pairwise overlap region in this case.
    triple_overlap_area = overlap_area
    
    # Total shadow area = Sum(A_i) - Sum(A_i n A_j) + (A1 n A2 n A3)
    total_shadow_area = sum_of_shadows - sum_of_pairwise_overlaps + triple_overlap_area

    # Step 5: Calculate the final area.
    final_area = initial_area - house_area + total_shadow_area

    print("--- Calculation Steps ---")
    print(f"1. Initial reachable area (no obstacles): {initial_area}")
    print(f"2. Area of the house (obstacle): {house_area}")
    print(f"3. Sum of individual shadow areas (ignoring overlaps): {sum_of_shadows}")
    print(f"4. Total shadow area (using Inclusion-Exclusion): {sum_of_shadows} - {sum_of_pairwise_overlaps} + {triple_overlap_area} = {total_shadow_area}")
    print("\n--- Final Result ---")
    print("The final reachable area is the Initial Area - House Area + Total Shadow Area.")
    print(f"Final Area = {initial_area} - {house_area} + {total_shadow_area} = {final_area}")

solve_horse_area()