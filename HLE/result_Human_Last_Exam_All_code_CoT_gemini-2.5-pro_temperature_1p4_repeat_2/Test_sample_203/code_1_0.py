import math

def solve_triangle_problem():
    """
    This function solves the geometry problem step-by-step and prints the result.
    """
    # Known values from the problem statement
    h_a = 12 # Height AD
    l_a = 13 # Angle bisector AE

    # --- Derivations explained in print statements ---

    print("Step 1: A coordinate system is established with D at (0,0) and A at (0,12).")

    # Step 2: Calculate the length of segment DE
    # In the right triangle ADE, AE^2 = AD^2 + DE^2
    DE_sq = l_a**2 - h_a**2
    DE = math.sqrt(DE_sq)
    print(f"Step 2: The length of segment DE is calculated from right triangle ADE: sqrt({l_a}^2 - {h_a}^2) = {DE}")

    # Step 3: A key relationship between the segments x=BD and y=CD is derived from angle properties.
    # This relationship is: 119*(y - x) = 1440 + 10*x*y
    print("Step 3: A key relationship is derived: 119*(y - x) = 1440 + 10*x*y")

    # Step 4: The condition for angle A to be acute (< 90 degrees) is xy < 144.
    # For a valid triangle, we also have x>0, y>0, so xy > 0.
    print("Step 4: The condition for angle A being acute gives the inequality: 0 < x*y < 144")

    # Step 5: The range for u = y - x is derived from the previous steps.
    # 0 < (119*u - 1440)/10 < 144  => 1440 < 119*u < 2880
    u_lower_bound = 1440 / 119
    u_upper_bound = 2880 / 119
    print(f"Step 5: This leads to the range for u = y - x: {u_lower_bound:.4f} < u < {u_upper_bound:.4f}")

    # Step 6: The median m is related to u by m^2 = 12^2 + (u/2)^2
    print("Step 6: The median m is related to u via m^2 = 12^2 + (u/2)^2")

    # Step 7: Calculate the final range for m
    m_sq_lower = h_a**2 + (u_lower_bound / 2)**2
    m_sq_upper = h_a**2 + (u_upper_bound / 2)**2

    m_lower = math.sqrt(m_sq_lower)
    m_upper = math.sqrt(m_sq_upper)
    
    # --- Final Answer ---

    print("\n--- Final Answer ---")
    print("The final equation describing the range of values for m is:")
    
    # Exact values for the upper bound
    m_upper_num = 12 * 169
    m_upper_den = 119

    # Exact values for the lower bound
    m_lower_num_part1 = 12
    m_lower_num_part2_sq = 119**2 + 60**2 # This is 17761
    m_lower_den = 119
    
    # Output the inequality with each number clearly stated as requested
    print(f"({m_lower_num_part1} * sqrt({m_lower_num_part2_sq}) / {m_lower_den}) < m < ({m_upper_num} / {m_upper_den})")
    print(f"As a decimal approximation:")
    print(f"{m_lower:.4f} < m < {m_upper:.4f}")
    
solve_triangle_problem()