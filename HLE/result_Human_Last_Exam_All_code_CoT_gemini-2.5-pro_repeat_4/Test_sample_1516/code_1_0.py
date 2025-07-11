import math

def solve_parliament_design():
    """
    Calculates the maximum integer value for K in the parliament design.
    """
    # Step 1: Determine the Seating Layout
    total_members = 791
    num_sections = 61
    # Each row in each section has one person
    num_rows_N = total_members // num_sections
    
    print(f"The number of rows per section is N = {num_rows_N}")

    # Step 2: Calculate Radial Distances
    r1 = 3.0  # Radial distance of the first row in meters
    row_depth = 1.5  # Radial distance between consecutive rows in meters

    # r_k = r1 + (k - 1) * d
    r2 = r1 + (2 - 1) * row_depth
    r_N = r1 + (num_rows_N - 1) * row_depth

    print(f"The radial distance of the first row is r_1 = {r1} m")
    print(f"The radial distance of the second row (the critical blocker) is r_2 = {r2} m")
    print(f"The radial distance of the last row (the viewer) is r_{num_rows_N} = {r_N} m")
    print("-" * 20)
    
    # Step 3, 4, 5: The derived visibility constraint is K < 2 * (r_2 - r_1) * (r_N + r_1)
    # This is the most restrictive case, considering a viewer in the last row (r_N),
    # a speaker in the opposite first row (-r_1), and a blocker in the opposite second row (-r_2).
    
    print("The derived visibility constraint, which must hold for the design to be stable, is:")
    print("K < 2 * (r_2 - r_1) * (r_N + r_1)")
    print()
    print("Substituting the calculated values into the equation:")
    
    # Step 6: Calculate the Maximum Integer K
    r2_minus_r1 = r2 - r1
    rN_plus_r1 = r_N + r1
    
    # Output each number in the final equation
    print(f"K < 2 * ({r2} - {r1}) * ({r_N} + {r1})")
    print(f"K < 2 * ({r2_minus_r1}) * ({rN_plus_r1})")
    
    upper_bound = 2 * r2_minus_r1 * rN_plus_r1
    
    print(f"K < {upper_bound}")
    
    # K must be an integer and strictly less than the upper bound.
    # The maximum integer K is floor(upper_bound - epsilon), which is upper_bound - 1
    # if upper_bound is an integer.
    max_K = math.ceil(upper_bound) - 1

    print("-" * 20)
    print(f"The maximum integer value K can take is {max_K}.")
    
    # Final answer format
    print(f"\n<<<{max_K}>>>")

solve_parliament_design()