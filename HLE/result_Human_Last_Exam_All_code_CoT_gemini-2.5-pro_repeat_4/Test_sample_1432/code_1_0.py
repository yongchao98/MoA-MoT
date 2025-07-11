def solve_chess_problem():
    """
    Calculates the minimum number of chess pieces to remove from a 7x8 board
    to prevent any sequence of 5 or more connected pieces in a straight line.
    """

    # Define the board dimensions.
    rows = 7
    cols = 8

    # The optimal strategy is to remove all pieces from one full row and one full column.
    # The number of pieces in a row is equal to the number of columns.
    pieces_in_a_row = cols
    
    # The number of pieces in a column is equal to the number of rows.
    pieces_in_a_column = rows

    # When removing a full row and a full column, the piece at their
    # intersection is counted twice, so we subtract 1.
    overlap = 1

    # Calculate the minimum number of pieces to be removed.
    min_removed = pieces_in_a_row + pieces_in_a_column - overlap

    # Print the explanation and the final equation with the result.
    print(f"The minimum number of pieces to be removed is calculated by removing one full row and one full column.")
    print(f"The equation is: (pieces in a row) + (pieces in a column) - (overlapping piece)")
    print(f"Result: {pieces_in_a_row} + {pieces_in_a_column} - {overlap} = {min_removed}")

solve_chess_problem()