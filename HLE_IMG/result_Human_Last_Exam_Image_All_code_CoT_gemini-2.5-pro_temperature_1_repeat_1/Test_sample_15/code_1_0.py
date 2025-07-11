def solve_puzzle():
    """
    Solves the visual puzzle by analyzing patterns in shapes and dots.
    """
    # Step 1: Analyze the shape pattern.
    # The shape is consistent per row. Row 1: Circle, Row 2: Square, Row 3: Triangle.
    # The missing piece is in Row 3, so its shape must be a Triangle.
    final_shape = "Triangle"
    print("Step 1: Analyzing the shapes.")
    print("The shape is consistent in each row.")
    print("Row 1 has Circles, Row 2 has Squares, Row 3 has Triangles.")
    print(f"Therefore, the missing shape is a {final_shape}.\n")

    # Step 2: Analyze the dot pattern by summing dots in each row.
    # Dots in the grid:
    # Row 1: [0, 4, 4]
    # Row 2: [0, 2, 2]
    # Row 3: [0, 2, ?]
    row1_dots = [0, 4, 4]
    row2_dots = [0, 2, 2]
    row3_dots_known = [0, 2]

    row1_sum = sum(row1_dots)
    row2_sum = sum(row2_dots)

    print("Step 2: Analyzing the dot pattern.")
    print(f"The sum of dots in Row 1 is: {row1_dots[0]} + {row1_dots[1]} + {row1_dots[2]} = {row1_sum}")
    print(f"The sum of dots in Row 2 is: {row2_dots[0]} + {row2_dots[1]} + {row2_dots[2]} = {row2_sum}")
    print("The total number of dots per row follows a geometric progression (halving each time).")
    
    # The pattern is 8 -> 4 -> ?
    # The ratio is 4 / 8 = 0.5
    expected_row3_sum = row2_sum / 2
    print(f"The expected sum for Row 3 is {row2_sum} / 2 = {int(expected_row3_sum)}.\n")

    # Step 3: Calculate the number of dots in the missing square.
    print("Step 3: Calculating the missing dot count.")
    # Expected Row 3 sum = known dots + missing dots
    # 2 = 0 + 2 + ?
    missing_dots = expected_row3_sum - sum(row3_dots_known)
    print(f"The equation for Row 3 is: {row3_dots_known[0]} + {row3_dots_known[1]} + ? = {int(expected_row3_sum)}")
    print(f"Solving for the missing value: ? = {int(expected_row3_sum)} - ({row3_dots_known[0]} + {row3_dots_known[1]}) = {int(missing_dots)}")
    
    print("\nConclusion: The missing image is a {} with {} dots.".format(final_shape, int(missing_dots)))

solve_puzzle()