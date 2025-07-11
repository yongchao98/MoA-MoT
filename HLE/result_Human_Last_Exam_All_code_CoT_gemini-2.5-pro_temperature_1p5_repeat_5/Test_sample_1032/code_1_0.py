import math

def solve():
    """
    Calculates the total reachable area for the horse.
    """
    # Length of the rope
    L = 7.0 / 2.0

    # Step 1: Calculate the area in Quadrants 1, 2, and 4.
    # This area is not affected by the house in Quadrant 3.
    # It consists of three triangles, each with area 1/2 * L^2.
    num_unobstructed_quadrants = 3
    area_one_quadrant_triangle = 0.5 * L**2
    area_q1_q2_q4 = num_unobstructed_quadrants * area_one_quadrant_triangle

    # Step 2: Calculate the extra area reachable in Quadrant 3.
    # The rope can pivot around house corners A(-2,0) and E(0,-2).

    # Distance to pivot A(-2,0)
    d_O_A = 2.0
    # Remaining rope length from pivot A
    rope_left_A = L - d_O_A
    
    # The area added by pivoting at A is a triangle in Q3.
    # Base of the triangle is from x=-3.5 to x=-0.5, so base = 3.0.
    # Height of the triangle is 1.5.
    base_A = 3.0
    height_A = 1.5
    area_from_A = 0.5 * base_A * height_A
    
    # Distance to pivot E(0,-2) is also 2.0, so the area added is identical.
    area_from_E = area_from_A
    
    # The two areas are disjoint, so the total area in Q3 is their sum.
    num_pivots_q3 = 2
    area_q3 = area_from_A + area_from_E
    
    # Step 3: Calculate the total area.
    total_area = area_q1_q2_q4 + area_q3

    # Print the breakdown of the final equation
    print("The final calculation is based on summing the areas from different regions.")
    print("Final Equation: (num_unobstructed_quadrants * area_one_quadrant) + (num_pivots_q3 * area_from_one_pivot)")
    print(f"Total Area = ({num_unobstructed_quadrants} * (1/2 * ({L})^2)) + ({num_pivots_q3} * (1/2 * {base_A} * {height_A}))")
    print(f"Total Area = {area_q1_q2_q4} + {area_q3} = {total_area}")

solve()
<<<22.875>>>