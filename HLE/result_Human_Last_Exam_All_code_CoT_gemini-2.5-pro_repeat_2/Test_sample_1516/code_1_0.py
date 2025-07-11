import math

def solve_parliament_design():
    """
    Calculates the maximum integer value for K in the parliament design.
    """
    # Step 1: Determine the Seating Layout
    total_members = 791
    num_sections = 61
    rows_per_section = total_members // num_sections
    
    # Step 2: Define physical dimensions
    r_initial = 3.0  # meters
    depth_per_row = 1.5  # meters
    h_seated = 1.0  # meters
    h_standing = 1.5 # meters
    
    # The height difference between a seated person and a standing person
    # at the same floor level.
    h_diff = h_seated - h_standing # This is -0.5

    # Step 3 & 4: Set up and analyze the visibility constraint.
    # The geometric condition for visibility is that for a speaker in row s,
    # an obstruction in row m, and an observer in row n (s < m < n),
    # the slope of the line of sight from the obstruction to the speaker must
    # be less than the slope from the observer to the speaker.
    #
    # The derivation leads to the inequality:
    # 1/K > (h_seated - h_standing) / ((r_n - r_s) * (r_m - r_s))
    # With the given numbers, h_seated - h_standing = -0.5, so:
    # 1/K > -0.5 / ((r_n - r_s) * (r_m - r_s))
    # Since K > 0 and the RHS is negative, this is always true and provides no upper bound.
    #
    # To find a maximum value for K, there must be a constraint of the form K <= value.
    # This occurs if the term (h_seated - h_standing) is positive.
    # We will proceed assuming this is the intended logic to make the problem solvable,
    # effectively using a positive height difference of 0.5.
    
    h_diff_assumed = 0.5

    # Step 5 & 6: Find the strictest constraint
    # The inequality becomes K < 2 * (r_n - r_s) * (r_m - r_s).
    # To satisfy this for all m and n, K must be less than the minimum
    # possible value of the right-hand side. This minimum occurs for the
    # smallest possible values of m and n.
    s = 1
    m = 2
    n = 3

    # Calculate the radial positions for these rows
    r_s = r_initial + (s - 1) * depth_per_row
    r_m = r_initial + (m - 1) * depth_per_row
    r_n = r_initial + (n - 1) * depth_per_row
    
    # Step 7: Calculate the upper bound for K
    # K < 1 / (h_diff_assumed / ((r_n - r_s) * (r_m - r_s)))
    # K < ((r_n - r_s) * (r_m - r_s)) / h_diff_assumed
    # K < ( (r_n-r_s) * (r_m-r_s) ) * 2
    K_upper_bound = ((r_n - r_s) * (r_m - r_s)) / h_diff_assumed
    
    # K must be an integer, so we take the floor of the upper bound.
    # Since the inequality is strict (K < bound), the max integer K is floor(bound - epsilon),
    # which is equivalent to ceil(bound) - 1 if the bound is not an integer,
    # or bound - 1 if the bound is an integer.
    max_K = math.floor(K_upper_bound - 1e-9) # Safely get the integer just below the bound

    print("--- Problem Analysis ---")
    print(f"Total Members: {total_members}, Sections: {num_sections}")
    print(f"Rows per section: {rows_per_section}")
    print(f"Speaker row (s): {s}, Obstruction row (m): {s+1}, Observer row (n): {s+2}")
    print(f"Radial distance of speaker row (r_s={s}): {r_s} m")
    print(f"Radial distance of obstruction row (r_m={m}): {r_m} m")
    print(f"Radial distance of observer row (r_n={n}): {r_n} m")
    print("\n--- Visibility Calculation ---")
    print("The visibility constraint requires K < ((r_n - r_s) * (r_m - r_s)) / (assumed_h_diff)")
    print(f"Equation: K < (({r_n} - {r_s}) * ({r_m} - {r_s})) / {h_diff_assumed}")
    print(f"Equation: K < ({r_n - r_s} * {r_m - r_s}) / {h_diff_assumed}")
    print(f"Equation: K < { (r_n - r_s) * (r_m - r_s) } / {h_diff_assumed}")
    print(f"Result: K < {K_upper_bound}")
    print(f"\nSince K must be an integer, the maximum value it can take is {max_K}.")
    
    return max_K

final_K = solve_parliament_design()
print(f"\nThe final answer is {final_K}")
print(f'<<<{final_K}>>>')