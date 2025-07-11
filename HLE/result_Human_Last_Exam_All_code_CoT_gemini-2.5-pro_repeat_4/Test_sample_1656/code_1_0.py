import collections

def solve_braid_index():
    """
    Calculates the braid index of a knot from its grid diagram representation.
    """
    # Grid number n
    n = 7

    # 'O' and 'X' positions from the problem description
    o_positions = [(1, 1), (2, 7), (3, 4), (4, 5), (5, 3), (6, 6), (7, 2)]
    x_positions = [(1, 2), (2, 6), (3, 3), (4, 1), (5, 7), (6, 5), (7, 4)]

    # Store the row number for each column's 'O' and 'X'
    # We use dictionaries for clarity, mapping column number to row number.
    o_rows = dict(o_positions)
    x_rows = dict(x_positions)

    braid_index = 0
    equation_terms = []

    print("Calculating the braid index by counting columns where the 'O' marker is above the 'X' marker.")
    print("-" * 70)

    # Iterate through each column from 1 to n
    for i in range(1, n + 1):
        o_row = o_rows[i]
        x_row = x_rows[i]

        # Check if the 'O' is in a higher row than the 'X'
        if o_row > x_row:
            is_greater = True
            equation_terms.append("1")
            braid_index += 1
            comparison_symbol = ">"
        else:
            is_greater = False
            equation_terms.append("0")
            comparison_symbol = "<="

        print(f"Column {i}: O is in row {o_row}, X is in row {x_row}. Comparison: {o_row} {comparison_symbol} {x_row}. Adding '{equation_terms[-1]}' to sum.")

    print("-" * 70)

    # Display the final calculation as a sum
    equation_str = " + ".join(equation_terms)
    print(f"The braid index is the result of the following sum:")
    print(f"{equation_str} = {braid_index}")
    
    print(f"\nThe braid index of the corresponding knot is {braid_index}.")

solve_braid_index()