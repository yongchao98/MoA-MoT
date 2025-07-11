import math

def solve_parliament_design():
    """
    Calculates the maximum integer value of K for the parliament design based on visibility constraints.
    """
    # Step 1: Define constants from the problem description
    num_members = 791
    num_sections = 61
    initial_radius_r1 = 3.0  # meters
    row_depth = 1.5  # meters
    
    # As reasoned in the plan, for the problem to have a maximum K, we must assume
    # the seated members are effectively higher than the standing speaker.
    # This could be due to a lowered podium for the speaker.
    # h_standing = 1.0
    # h_seated = 1.5
    # The difference in height is what matters. h_seated - h_standing = 0.5m
    height_diff = 0.5 # Corresponds to (h_seated) - (h_standing) or -(h_standing - h_seated)

    # Step 2: Calculate the number of rows
    num_rows = num_members / num_sections
    print(f"Number of rows per section = {num_members} / {num_sections} = {int(num_rows)}")
    
    # Step 3: Identify the critical case for visibility
    # The constraint is most restrictive for the closest possible observer and blocker.
    # Observer in row i=3, Blocker in row j=2, Speaker in row 1.
    i = 3
    j = 2
    
    print(f"\nThe visibility constraint is most restrictive for an observer in row {i} and a blocker in row {j}.")
    
    # Step 4: Calculate the radial positions for the critical case
    # r_i = r_1 + (i-1) * depth
    r1 = initial_radius_r1
    r2 = r1 + (j - 1) * row_depth
    r3 = r1 + (i - 1) * row_depth
    
    print(f"The radius of the first row (r1) = {r1} m")
    print(f"The radius of the second row (r2) = {r1} + (2-1)*{row_depth} = {r2} m")
    print(f"The radius of the third row (r3) = {r1} + (3-1)*{row_depth} = {r3} m")

    # Step 5: Derive and calculate the upper bound for K
    # The derived inequality for K is: K < 2 * (r_i - r1) * (r_j - r1) / height_diff
    # Here height_diff is 0.5 as it comes from a term like `... > height_diff * ...`
    # Our algebra based on the height assumption led to: K < 2 * (r_i - r1) * (r_j - r1)
    # The factor of 2 in the expression is effectively 1 / height_diff, so the formulas match.
    
    r3_minus_r1 = r3 - r1
    r2_minus_r1 = r2 - r1
    
    limit_k = 2 * r3_minus_r1 * r2_minus_r1
    
    print("\nThe inequality derived from the visibility condition is:")
    print(f"K < 2 * (r{i} - r1) * (r{j} - r1)")
    print(f"K < 2 * ({r3} - {r1}) * ({r2} - {r1})")
    print(f"K < 2 * {r3_minus_r1} * {r2_minus_r1}")
    print(f"K < {limit_k}")
    
    # Step 6: Find the maximum integer K that satisfies this condition
    # Since K must be strictly less than the limit, the max integer K is floor(limit - epsilon).
    max_k = math.floor(limit_k - 1e-9)
    
    print(f"\nSince K must be an integer, the maximum value it can take is {max_k}.")

solve_parliament_design()
<<<8>>>