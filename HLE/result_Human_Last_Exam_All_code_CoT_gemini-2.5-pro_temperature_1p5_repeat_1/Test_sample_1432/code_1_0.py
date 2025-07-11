def solve_chess_puzzle():
    """
    Calculates the minimum number of pieces to remove from a 7x8 board
    to prevent any 5-in-a-row sequences.
    
    The strategy is to remove all pieces from a single row and a single column.
    """
    
    rows = 7
    cols = 8
    
    # The number of pieces in a single column is equal to the number of rows.
    pieces_in_one_col = rows
    
    # The number of pieces in a single row is equal to the number of columns.
    pieces_in_one_row = cols
    
    # When we remove one row and one column, the piece at their intersection
    # is counted in both removals, so we subtract 1.
    overlap = 1
    
    min_removed = pieces_in_one_row + pieces_in_one_col - overlap
    
    print("The strategy is to remove all pieces from one full row and one full column.")
    print("This guarantees that no sequence of 5 pieces can be formed in any direction.")
    print("The calculation is based on the board dimensions (7 rows, 8 columns):")
    print(f"Pieces in a row ({pieces_in_one_row}) + Pieces in a column ({pieces_in_one_col}) - Overlapping piece ({overlap}) = Total removed")
    print(f"{pieces_in_one_row} + {pieces_in_one_col} - {overlap} = {min_removed}")

solve_chess_puzzle()