def is_square_attacked(square, white_king, white_queen, white_pawn):
    """
    Checks if a given square is attacked by any of the white pieces.
    Note: For this specific problem, we assume no pieces are blocking the queen's path.
    """
    s_col, s_row = square
    wk_col, wk_row = white_king
    wq_col, wq_row = white_queen
    wp_col, wp_row = white_pawn

    # 1. Check for attack by the White King
    if abs(s_col - wk_col) <= 1 and abs(s_row - wk_row) <= 1:
        # The king cannot attack its own square, but this is handled by the main loop
        return True

    # 2. Check for attack by the White Queen (Rook + Bishop moves)
    # Rank and file (Rook)
    if s_col == wq_col or s_row == wq_row:
        return True
    # Diagonals (Bishop)
    if abs(s_col - wq_col) == abs(s_row - wq_row):
        return True

    # 3. Check for attack by the White Pawn
    # A white pawn at (col, row) attacks (col-1, row+1) and (col+1, row+1)
    if s_row == wp_row + 1 and abs(s_col - wp_col) == 1:
        return True

    return False

def solve_chess_puzzle():
    """
    Verifies the solution to the chess stalemate problem and prints the result.
    """
    # Piece positions are 0-indexed: (column, row)
    # 'a' to 'h' -> 0 to 7
    # '1' to '8' -> 0 to 7
    white_king_pos = (6, 4)  # g5
    white_queen_pos = (3, 3) # d4
    white_pawn_pos = (4, 1)  # e2
    black_king_pos = (7, 0)  # h1

    white_pieces_pos = [white_king_pos, white_queen_pos, white_pawn_pos]
    
    unattacked_squares = []

    for col in range(8):
        for row in range(8):
            square_to_check = (col, row)
            
            # A piece does not attack the square it occupies
            if square_to_check in white_pieces_pos:
                continue

            if not is_square_attacked(square_to_check, white_king_pos, white_queen_pos, white_pawn_pos):
                # Convert 0-indexed coordinates to chess notation for printing
                col_char = chr(ord('a') + col)
                row_char = str(row + 1)
                unattacked_squares.append(f"{col_char}{row_char}")

    print("Verifying the position: White King at g5, Queen at d4, Pawn at e2. Black King at h1.")
    print(f"Squares that are not attacked by any white piece: {unattacked_squares}")
    
    if len(unattacked_squares) == 1 and unattacked_squares[0] == "h1":
        print("\nVerification successful. Only h1 is unattacked, resulting in a stalemate for the Black King.")
        queen_points = 9
        pawn_points = 1
        total_points = queen_points + pawn_points
        print("\nThe smallest number of points is calculated as follows:")
        print(f"Queen ({queen_points}) + Pawn ({pawn_points}) = {total_points}")
    else:
        print("\nThe provided position is not a valid solution.")

solve_chess_puzzle()
<<<10>>>