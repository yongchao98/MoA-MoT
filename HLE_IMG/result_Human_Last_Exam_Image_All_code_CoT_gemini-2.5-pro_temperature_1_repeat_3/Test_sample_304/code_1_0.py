def solve_group_exponent_sum():
    """
    This script calculates the sum of group exponents for each column of the provided 4x4 grid of graphs.
    """

    # Exponents of the identified groups for each graph V_1 to V_16
    # The grid is represented as a list of lists (rows of columns)
    exponent_grid = [
        # Row 1: V1, V2, V3, V4
        [30, 12, 6, 4],
        # Row 2: V5, V6, V7, V8
        [8, 2, 6, 3],
        # Row 3: V9, V10, V11, V12
        [30, 3, 8, 6],
        # Row 4: V13, V14, V15, V16
        [4, 12, 2, 6]
    ]

    # Initialize column sums
    s1, s2, s3, s4 = 0, 0, 0, 0

    # Extract exponents for each column
    col1_exponents = [row[0] for row in exponent_grid]
    col2_exponents = [row[1] for row in exponent_grid]
    col3_exponents = [row[2] for row in exponent_grid]
    col4_exponents = [row[3] for row in exponent_grid]

    # Calculate sums
    s1 = sum(col1_exponents)
    s2 = sum(col2_exponents)
    s3 = sum(col3_exponents)
    s4 = sum(col4_exponents)

    # Print the breakdown of each sum
    print("Column Sum Calculations:")
    print(f"S1 = {col1_exponents[0]} + {col1_exponents[1]} + {col1_exponents[2]} + {col1_exponents[3]} = {s1}")
    print(f"S2 = {col2_exponents[0]} + {col2_exponents[1]} + {col2_exponents[2]} + {col2_exponents[3]} = {s2}")
    print(f"S3 = {col3_exponents[0]} + {col3_exponents[1]} + {col3_exponents[2]} + {col3_exponents[3]} = {s3}")
    print(f"S4 = {col4_exponents[0]} + {col4_exponents[1]} + {col4_exponents[2]} + {col4_exponents[3]} = {s4}")

    # Print the final result in the specified format
    final_answer = f"{{{s1}, {s2}, {s3}, {s4}}}"
    print("\nFinal ordered list of column sums:")
    print(final_answer)

# Execute the function
solve_group_exponent_sum()