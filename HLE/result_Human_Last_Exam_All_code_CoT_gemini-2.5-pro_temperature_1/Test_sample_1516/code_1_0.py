import math

def solve_parliament_design():
    """
    Calculates the maximum integer value for K in the parliament design.
    """
    # Step 1: Define initial parameters from the problem description
    total_members = 791
    num_sections = 61
    person_depth = 1.5  # meters
    initial_radius = 3.0  # meters (r_0)
    
    # Calculate the number of rows per section
    # Since each row seats one person, this is total members / sections.
    num_rows = total_members // num_sections
    
    print(f"Number of members: {total_members}")
    print(f"Number of sections: {num_sections}")
    print(f"Each section will have {num_rows} rows.\n")

    # Step 2: Formulate the visibility constraint
    # The parliament floor follows the equation: h = r^2 / K
    # A person speaking stands in the first row (row 0), at radius r_0.
    # The speaker's feet are at floor height: h_feet(0) = r_0^2 / K
    # A person sits in row 'i'. Their head is at height H_head(i) = (r_i^2 / K) + 1.0
    #
    # To ensure full visibility, the line of sight from an observer in row 'i+1' to the
    # speaker's feet must clear the head of the blocker in row 'i'.
    # This leads to the inequality: K <= (r_i - r_0) * (r_{i+1} - r_0)
    
    print("The visibility constraint requires that the line of sight from any member to the speaker")
    print("is not obstructed. We use the most restrictive model: the line of sight from an")
    print("observer's eyes to the speaker's feet must clear the head of any person in between.\n")

    print("This leads to the following inequality for the constant K:")
    print("K <= (r_i - r_0) * (r_{i+1} - r_0)\n")

    # Step 3: Find the most restrictive case for the inequality
    # The term on the right side, (r_i - r_0)*(r_{i+1} - r_0), grows as 'i' increases.
    # Therefore, the most restrictive condition (lowest upper bound for K) is for the first few rows, at i=1.
    # We must satisfy: K <= (r_1 - r_0) * (r_2 - r_0)
    
    print("This inequality must hold for all rows. The most restrictive case is for the rows")
    print("closest to the front (i=1), which provides the lowest upper bound for K.")
    print("So, we evaluate the inequality for i=1:\n")
    print("K <= (r_1 - r_0) * (r_2 - r_0)\n")

    # Step 4: Calculate the values and solve for K
    r_0 = initial_radius
    r_1 = initial_radius + 1 * person_depth
    r_2 = initial_radius + 2 * person_depth
    
    delta_r1 = r_1 - r_0
    delta_r2 = r_2 - r_0
    
    print("Calculating the radial distances:")
    print(f"r_0 (radius of row 0) = {r_0} m")
    print(f"r_1 (radius of row 1) = {r_1} m")
    print(f"r_2 (radius of row 2) = {r_2} m\n")

    print("Plugging these values into the inequality:")
    print(f"K <= ({r_1} - {r_0}) * ({r_2} - {r_0})")
    print(f"K <= {delta_r1} * {delta_r2}")
    
    max_K_float = delta_r1 * delta_r2
    print(f"K <= {max_K_float}\n")
    
    # Step 5: Determine the maximum integer value for K
    max_K_int = math.floor(max_K_float)
    
    print(f"Since K must be an integer, the maximum possible value for K is {max_K_int}.")
    
    print(f"\n<<<4>>>")

solve_parliament_design()