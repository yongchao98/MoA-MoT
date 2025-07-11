import math

def solve_parliament_design():
    """
    Calculates the maximum integer value for the parameter K in the parliament design.
    """
    # Step 1: Define the geometric parameters
    num_members = 791
    num_sections = 61
    num_rows = num_members // num_sections
    
    r_initial = 3.0  # meters
    row_depth = 1.5  # meters
    
    # Calculate radius of the first and last rows
    r1 = r_initial
    # r_i = r_initial + (i - 1) * row_depth
    r13 = r_initial + (num_rows - 1) * row_depth
    
    print(f"Number of rows per section: {num_rows}")
    print(f"Radius of the first row (r1): {r1} m")
    print(f"Radius of the last row (r13): {r13} m\n")

    # Step 2 & 3: Formulate and explain the critical constraint.
    # The height of a seated person in row 'i' is Hi = ri^2/K + 1.
    # The height of the standing speaker in the first row is Hs = r1^2/K + 1.5.
    # The problem implies a constraint beyond simple line-of-sight visibility, related to parliamentary ethos and status.
    # A flat hierarchy is desired, so we can infer that members in the back should not feel they are looking down on the speaker.
    # This leads to the constraint that the head height of a member in the last row must be at least the head height of the speaker.
    # H_seated_13 >= H_speaker
    # r13^2/K + 1 >= r1^2/K + 1.5
    
    # Step 4: Solve the inequality for K
    # r13^2/K - r1^2/K >= 1.5 - 1
    # (r13^2 - r1^2) / K >= 0.5
    # Since K must be positive (paraboloid opening upwards), we can multiply by K.
    # r13^2 - r1^2 >= 0.5 * K
    # K <= (r13^2 - r1^2) / 0.5
    # K <= 2 * (r13^2 - r1^2)
    
    r1_sq = r1**2
    r13_sq = r13**2
    
    max_K_float = 2 * (r13_sq - r1_sq)
    
    # K must be an integer
    max_K_integer = math.floor(max_K_float)

    print("The derived inequality based on the 'status' constraint is:")
    print(f"K <= 2 * (r13^2 - r1^2)")
    print(f"K <= 2 * ({r13_sq} - {r1_sq})")
    print(f"K <= 2 * ({r13_sq - r1_sq})")
    print(f"K <= {max_K_float}\n")
    print(f"Since K must be an integer, the maximum value K can take is {max_K_integer}.")
    
    # The final answer format
    # print(f"\n<<<{max_K_integer}>>>")

solve_parliament_design()