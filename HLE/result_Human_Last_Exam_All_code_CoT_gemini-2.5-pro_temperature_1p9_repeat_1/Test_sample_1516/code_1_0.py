import math

def solve_parliament_design():
    """
    Calculates the maximum integer value for K based on the visibility constraints
    for a paraboloid-shaped parliament hall.
    """

    # Step 1 & 2: Define the parameters from the problem description
    num_members = 791
    num_sections = 61
    num_rows = num_members // num_sections
    
    initial_radius = 3.0  # r_1 in meters
    depth_per_row = 1.5   # in meters
    
    print("Parliament Design Analysis")
    print("==========================")
    print(f"Number of members: {num_members}")
    print(f"Number of sections: {num_sections}")
    print(f"Number of rows per section: {num_rows}\n")

    # Step 3 & 4: Identify the critical case for the visibility constraint.
    # The constraint is K < (r_n - r_1) * (r_m - r_1).
    # This constraint is tightest (the right side is smallest) for the smallest
    # possible values of n and m where an obstruction can occur.
    # This is when an observer in row n=3 is viewing the speaker in row 1,
    # with a potential obstruction from the person in row m=2.
    
    n_observer = 3
    m_obstruction = 2
    
    print("Identifying the Critical Case for Visibility:")
    print("The visibility constraint is derived from the condition that an observer's line of sight")
    print("to the speaker's feet is not blocked by anyone's head in an intermediate row.")
    print("This leads to the inequality: K < (r_n - r_1) * (r_m - r_1)")
    print(f"The most restrictive case is for the closest observer, n = {n_observer}, and closest obstruction, m = {m_obstruction}.\n")

    # Step 5: Calculate the radii for the critical rows.
    # r_n = initial_radius + (n - 1) * depth_per_row
    r1 = initial_radius
    r2 = initial_radius + (m_obstruction - 1) * depth_per_row
    r3 = initial_radius + (n_observer - 1) * depth_per_row
    
    print("Calculating Radii for Critical Rows:")
    print(f"Radius of speaker's row (r_1): {r1} m")
    print(f"Radius of obstruction row (r_2): {r2} m")
    print(f"Radius of observer's row (r_3): {r3} m\n")

    # Step 6: Apply the values to the inequality to find the bound for K.
    r_n_minus_r_1 = r3 - r1
    r_m_minus_r_1 = r2 - r1
    
    bound_for_K = r_n_minus_r_1 * r_m_minus_r_1
    
    print("Calculating the Upper Bound for K:")
    print(f"The inequality is: K < (r_3 - r_1) * (r_2 - r_1)")
    # We output each number in the final equation as requested
    print(f"Substituting values: K < ({r_n_minus_r_1}) * ({r_m_minus_r_1})")
    print(f"This simplifies to: K < {bound_for_K}\n")
    
    # The problem requires K to be an integer. Since K must be strictly less
    # than the bound, the maximum integer value is floor(bound - epsilon).
    # If the bound is an integer, K_max = bound - 1. Otherwise, K_max = floor(bound).
    if bound_for_K == int(bound_for_K):
        max_K = int(bound_for_K) - 1
    else:
        max_K = math.floor(bound_for_K)
        
    print("Determining the Final Answer:")
    print(f"Since K must be an integer and strictly less than {bound_for_K}, the maximum possible integer value for K is {max_K}.")
    
    return max_K

# Execute the function and print the final answer in the required format.
final_answer = solve_parliament_design()
print(f"\n<<<4>>>")
