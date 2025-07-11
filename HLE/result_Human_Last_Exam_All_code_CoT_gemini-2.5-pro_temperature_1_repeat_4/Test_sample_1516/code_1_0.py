import math

def solve_parliament_design():
    """
    Calculates the maximum integer value for K in the parliament design problem.
    """
    # Step 1: Define constants from the problem description.
    total_members = 791
    num_sections = 61
    initial_radius = 3.0  # meters
    row_depth = 1.5       # meters
    seated_height = 1.0   # meters
    standing_height = 1.5 # meters

    # Step 2: Calculate the number of rows.
    # Each row in each section has one person.
    num_rows = total_members // num_sections
    
    print(f"Number of members: {total_members}")
    print(f"Number of sections: {num_sections}")
    print(f"Calculated number of rows per section: {num_rows}\n")

    # Step 3: Define the visibility condition and derive the inequality for K.
    # The geometric condition for visibility for an observer in row n to see the speaker
    # in row 1 over the blocker in row n-1 simplifies to:
    # K < 2 * (r_{n-1} - r_1) * (r_n - r_1)
    # This must hold for n = 3, 4, ..., num_rows.
    # We need to find the most restrictive (minimum) limit on K.

    print("Deriving the constraint on K for each row...")
    
    min_k_limit = float('inf')

    # The first row's radial distance is fixed.
    r1 = initial_radius

    # Step 4: Iterate from n=3 to num_rows to find the tightest constraint.
    # The check starts at n=3 because that's the first row with both an observer and a blocker.
    for n in range(3, num_rows + 1):
        # Calculate the radial distance for the relevant rows
        r_n_minus_1 = initial_radius + (n - 2) * row_depth
        r_n = initial_radius + (n - 1) * row_depth

        # Calculate the upper limit for K for this specific configuration
        k_limit = 2 * (r_n_minus_1 - r1) * (r_n - r1)
        
        print(f"For an observer in row n={n}:")
        print(f"  - Blocker is in row {n-1} at radius r_{n-1} = {r_n_minus_1:.1f}m")
        print(f"  - Observer is in row {n} at radius r_{n} = {r_n:.1f}m")
        print(f"  - The condition is K < {k_limit:.1f}")

        # Update the minimum limit found so far
        if k_limit < min_k_limit:
            min_k_limit = k_limit

    print(f"\nThe most restrictive constraint found is K < {min_k_limit:.1f}\n")

    # Step 5: Determine the maximum integer value for K.
    # Since K must be an integer and K must be strictly less than the minimum limit,
    # we take the floor of the limit and subtract 1 if the limit is an integer,
    # or just take the floor if it's not.
    max_k = math.floor(min_k_limit - 0.000001)

    print(f"As K must be an integer, the maximum possible value for K is {max_k}.")
    print("\nThe final equation for the parliament's structure is:")
    print("h", "=", "r^2", "/", max_k)

solve_parliament_design()
<<<8>>>