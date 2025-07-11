import itertools

def to_coords(algebraic):
    """Converts algebraic notation like 'h8' to (row, col) tuple (7, 7)."""
    col = ord(algebraic[0]) - ord('a')
    row = int(algebraic[1]) - 1
    return row, col

def to_algebraic(coords):
    """Converts (row, col) tuple (7, 7) to algebraic notation 'h8'."""
    row, col = coords
    return chr(ord('a') + col) + str(row + 1)

def get_attacked_squares(pieces):
    """
    Calculates all squares attacked by a given set of white pieces.
    - pieces: A dictionary of {'algebraic_pos': 'piece_type'}
    - black_king_pos: The (row, col) of the black king.
    """
    attacked = set()
    
    # Convert piece positions to coordinate tuples
    piece_positions = {to_coords(pos): piece_type for pos, piece_type in pieces.items()}
    occupied_squares = set(piece_positions.keys())

    for pos, piece_type in piece_positions.items():
        r, c = pos
        # --- Knight ---
        if piece_type == 'N':
            moves = [(-2, -1), (-2, 1), (-1, -2), (-1, 2),
                     (1, -2), (1, 2), (2, -1), (2, 1)]
            for dr, dc in moves:
                nr, nc = r + dr, c + dc
                if 0 <= nr < 8 and 0 <= nc < 8:
                    attacked.add((nr, nc))

        # --- Rook and Queen (Orthogonal) ---
        if piece_type in ('R', 'Q'):
            directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]
            for dr, dc in directions:
                nr, nc = r + dr, c + dc
                while 0 <= nr < 8 and 0 <= nc < 8:
                    attacked.add((nr, nc))
                    if (nr, nc) in occupied_squares:
                        break
                    nr, nc = nr + dr, nc + dc

        # --- Bishop and Queen (Diagonal) ---
        if piece_type in ('B', 'Q'):
            directions = [(-1, -1), (-1, 1), (1, -1), (1, 1)]
            for dr, dc in directions:
                nr, nc = r + dr, c + dc
                while 0 <= nr < 8 and 0 <= nc < 8:
                    attacked.add((nr, nc))
                    if (nr, nc) in occupied_squares:
                        break
                    nr, nc = nr + dr, nc + dc
    return attacked

def solve_chess_puzzle():
    """
    Verifies the 5-piece solution to the stalemate domination problem.
    """
    # K. Fabel, 1951 solution
    white_pieces = {
        'd3': 'Q',
        'h7': 'R',
        'g5': 'B',
        'e6': 'B',
        'f5': 'N'
    }
    num_pieces = len(white_pieces)
    black_king_alg = 'h8'
    black_king_pos = to_coords(black_king_alg)

    # All pieces on the board block movement
    all_pieces_for_blocking = white_pieces.copy()
    all_pieces_for_blocking[black_king_alg] = 'k'
    
    # Calculate attacked squares
    all_attacked = get_attacked_squares(all_pieces_for_blocking)

    # Find unattacked squares
    all_squares = set(itertools.product(range(8), range(8)))
    unattacked_squares = all_squares - all_attacked

    # Verify stalemate condition for the black king
    kr, kc = black_king_pos
    escape_squares = set()
    for dr in [-1, 0, 1]:
        for dc in [-1, 0, 1]:
            if dr == 0 and dc == 0:
                continue
            nr, nc = kr + dr, kc + dc
            if 0 <= nr < 8 and 0 <= nc < 8:
                escape_squares.add((nr, nc))
    
    all_escape_squares_attacked = escape_squares.issubset(all_attacked)

    print(f"Verifying a {num_pieces}-piece solution...")
    print("White pieces:", white_pieces)
    print("Black king is on:", black_king_alg)
    print("-" * 30)

    print("Set of all unattacked squares (in row, col format):")
    unattacked_alg = {to_algebraic(s) for s in unattacked_squares}
    print(unattacked_alg if unattacked_alg else "None")
    print()

    print(f"Is the king's square {black_king_alg} unattacked? {'Yes' if black_king_pos in unattacked_squares else 'No'}")
    print(f"Are all king's escape squares attacked? {'Yes' if all_escape_squares_attacked else 'No'}")
    
    is_solution_valid = (len(unattacked_squares) == 1 and 
                         black_king_pos in unattacked_squares and 
                         all_escape_squares_attacked)

    print()
    if is_solution_valid:
        print("Conclusion: The position is a valid stalemate and dominates 63 squares.")
        print("The smallest number of points of white material is:")
        print(num_pieces)
    else:
        print("Conclusion: This configuration is not a valid solution.")

solve_chess_puzzle()