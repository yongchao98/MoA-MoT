import math

def solve_parliament_design():
    """
    Calculates the maximum integer value for K based on visibility constraints.
    """
    # Step 1 & 2: Define the parliament layout parameters
    initial_radius_r1 = 3.0  # meters
    row_depth = 1.5          # meters
    
    # The radius of the second row (the first potential obstruction)
    radius_r2 = initial_radius_r1 + row_depth
    
    # Step 6 & 7: Apply the visibility constraint to find the limit for K.
    # The derived condition is K < (r - r1)^2.
    # To satisfy this for all rows, K must be less than the minimum value of this expression.
    # The minimum occurs at the smallest possible radius r, which is r2.
    # u represents the distance from the first row, (r - r1)
    min_u = radius_r2 - initial_radius_r1
    k_limit = min_u**2
    
    # Step 8: The maximum integer value for K is the floor of this limit.
    max_k = math.floor(k_limit)
    # The problem implies K < limit, so if limit is an integer, K must be limit-1.
    # For example, if K < 9.0, max integer K is 8. If K < 2.25, max integer K is 2.
    if k_limit == float(max_k):
        max_k = max_k -1
    else:
        max_k = int(max_k)

    print("The visibility constraint is derived as: K < (r - r1)^2")
    print("To hold for all rows, K must be smaller than the minimum possible value of (r - r1)^2.")
    print(f"The minimum value occurs for the second row (r = r2).")
    print(f"r1 = {initial_radius_r1}")
    print(f"r2 = {radius_r2}")
    print(f"The calculation is: K < ({radius_r2} - {initial_radius_r1})^2")
    print(f"K < {min_u}^2")
    print(f"K < {k_limit}")
    print(f"\nThe maximum integer value for K is therefore {max_k}.")

solve_parliament_design()
<<<2>>>