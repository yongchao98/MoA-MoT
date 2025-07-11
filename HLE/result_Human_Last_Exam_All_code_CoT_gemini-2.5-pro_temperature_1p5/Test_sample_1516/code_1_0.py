import math

def solve_parliament_design():
    """
    Calculates the maximum integer value for the parameter K in the new parliament design.
    """

    # Step 1: Deconstruct the Parliament's Layout
    total_members = 791
    num_sections = 61
    rows_per_section = total_members / num_sections

    row_depth = 1.5  # in meters
    initial_radius = 3.0  # r_1, in meters

    # The problem implies there are 13 rows (791/61=13)
    # We only need the radius of the first two rows to find the tightest constraint.
    r1 = initial_radius
    r2 = initial_radius + (2 - 1) * row_depth

    print("Step 1: Parliament Layout")
    print(f"Number of rows per section: {total_members} / {num_sections} = {int(rows_per_section)}")
    print(f"Radius of the 1st row (r1): {r1} m")
    print(f"Radius of the 2nd row (r2): {r2} m")
    print("-" * 20)

    # Step 2 & 3: Model Heights and Establish Visibility Constraint
    # The floor height 'h' at radius 'r' is h = r^2 / K.
    # Speaker (row 1, standing) head height: H_1 = r1^2/K + 1.5
    # Seated member (row j > 1) head height: H_j = rj^2/K + 1.0
    #
    # The visibility constraint is interpreted as every seated member looking horizontally
    # or upwards at the speaker. This avoids anyone feeling 'higher' than the speaker.
    # Mathematically: Head_Height(seated_member) >= Head_Height(speaker)
    # H_j >= H_1  for all j > 1
    # rj^2/K + 1.0 >= r1^2/K + 1.5
    # (rj^2 - r1^2) / K >= 0.5
    # K <= 2 * (rj^2 - r1^2)
    #
    # To satisfy this for all rows, K must be less than or equal to the minimum value
    # of the expression. The minimum occurs at the smallest j, which is j=2.
    
    print("Step 2 & 3: Visibility Constraint")
    print("The key constraint is that every member must look horizontally or upwards at the speaker.")
    print("This means the head height of any seated member must be >= the speaker's head height.")
    print("This leads to the inequality: K <= 2 * (r_j^2 - r_1^2)")
    print("The tightest constraint occurs for the row closest to the speaker, which is row 2 (j=2).")
    print("-" * 20)

    # Step 4: Formulate and Solve the Inequality for the tightest constraint (j=2)
    # K <= 2 * (r2^2 - r1^2)
    
    r1_sq = r1**2
    r2_sq = r2**2
    
    max_k_float = 2 * (r2_sq - r1_sq)

    print("Step 4: Solving the Inequality")
    print(f"Using the tightest constraint with r1 = {r1} and r2 = {r2}:")
    print(f"K <= 2 * ({r2}^2 - {r1}^2)")
    print(f"K <= 2 * ({r2_sq} - {r1_sq})")
    print(f"K <= 2 * ({r2_sq - r1_sq})")
    print(f"K <= {max_k_float}")
    print("-" * 20)

    # Step 5: Determine the Maximum Integer Value for K
    max_k_int = math.floor(max_k_float)
    
    print("Step 5: Final Calculation")
    print("Since K must be an integer, we take the floor of the result.")
    print(f"The maximum integer value for K is {max_k_int}.")
    
    # Final Answer
    print("\nFinal Equation with values:")
    print(f"Max K = floor(2 * ({r2}^2 - {r1}^2))")
    print(f"Max K = floor({max_k_float})")
    print(f"Final Answer: {max_k_int}")


solve_parliament_design()
print("<<<22>>>")