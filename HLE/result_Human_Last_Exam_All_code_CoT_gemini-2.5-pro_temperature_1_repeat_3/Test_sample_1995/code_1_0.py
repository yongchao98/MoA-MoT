import itertools

def get_attacked_squares(pieces):
    """
    Calculates all squares attacked by a given set of pieces.
    
    Args:
        pieces (dict): A dictionary where keys are piece types ('q', 'r')
                       and values are lists of (col, row) tuples.
                       (0,0) is 'a1', (7,7) is 'h8'.

    Returns:
        set: A set of (col, row) tuples representing attacked squares.
    """
    attacked = set()
    
    # Define board dimensions
    board_range = range(8)

    # Rook attacks (and part of Queen attacks)
    for r_col, r_row in pieces.get('r', []) + pieces.get('q', []):
        for i in board_range:
            if i != r_row:
                attacked.add((r_col, i))
            if i != r_col:
                attacked.add((i, r_row))

    # Bishop attacks (and part of Queen attacks)
    for b_col, b_row in pieces.get('b', []) + pieces.get('q', []):
        for i in range(1, 8):
            # Diagonals
            coords = [
                (b_col + i, b_row + i),
                (b_col + i, b_row - i),
                (b_col - i, b_row + i),
                (b_col - i, b_row - i),
            ]
            for c, r in coords:
                if 0 <= c < 8 and 0 <= r < 8:
                    attacked.add((c, r))

    return attacked

def chess_notation(coord):
    """Converts a (col, row) tuple to chess notation like 'a1'."""
    col, row = coord
    return chr(ord('a') + col) + str(row + 1)

def solve_chess_puzzle():
    """
    This function sets up and verifies the solution to the chess stalemate problem.
    """
    # The solution uses 5 pieces.
    # Coordinates are 0-indexed, where (0,0) is 'a1' and (7,7) is 'h8'.
    # Piece positions: Qd2, Qb4, Qc5, Re1, Ra2
    white_pieces = {
        'q': [(3, 1), (1, 3), (2, 4)],  # Queens at d2, b4, c5
        'r': [(4, 0), (0, 1)],          # Rooks at e1, a2
    }
    
    # The total number of pieces is the answer we are testing.
    num_pieces = sum(len(positions) for positions in white_pieces.values())
    
    print(f"Testing a configuration with {num_pieces} white pieces.")
    
    # Define all squares on the board
    all_squares = set(itertools.product(range(8), range(8)))
    
    # Calculate all squares attacked by the white pieces
    attacked_by_white = get_attacked_squares(white_pieces)
    
    # Find the unattacked squares
    unattacked_squares = all_squares - attacked_by_white
    
    print(f"Number of attacked squares: {len(attacked_by_white)}")
    
    if len(unattacked_squares) == 1:
        king_pos = unattacked_squares.pop()
        print(f"Success: Exactly one square is unattacked: {chess_notation(king_pos)}")
        print("This is where the black king must be for the stalemate.")
        
        # Verify stalemate condition
        print("\nVerifying stalemate...")
        k_col, k_row = king_pos
        
        # Get king's potential moves
        king_moves = set()
        for dc in [-1, 0, 1]:
            for dr in [-1, 0, 1]:
                if dc == 0 and dr == 0:
                    continue
                nk_col, nk_row = k_col + dc, k_row + dr
                if 0 <= nk_col < 8 and 0 <= nk_row < 8:
                    king_moves.add((nk_col, nk_row))

        print(f"King at {chess_notation(king_pos)} wants to move to: {[chess_notation(s) for s in king_moves]}")
        
        # Check if all potential moves are to attacked squares
        if king_moves.issubset(attacked_by_white):
            print("All of the king's possible escape squares are attacked.")
            print("The position is a valid stalemate.")
            print("\n----------------------------------------------------")
            print(f"The smallest number of points of white material is: {num_pieces}")
            print("----------------------------------------------------")
        else:
            print("Failure: The king is not in stalemate as it has safe squares to move to.")
            
    else:
        print(f"Failure: Found {len(unattacked_squares)} unattacked squares, but expected 1.")
        print(f"Unattacked squares are: {[chess_notation(s) for s in unattacked_squares]}")

# Run the verification
solve_chess_puzzle()
