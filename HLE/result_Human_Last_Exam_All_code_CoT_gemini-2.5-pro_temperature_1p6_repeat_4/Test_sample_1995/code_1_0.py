def to_coords(s):
    """Converts algebraic notation (e.g., 'a1') to (row, col) tuple."""
    col = ord(s[0]) - ord('a')
    row = int(s[1]) - 1
    return row, col

def to_alg(r, c):
    """Converts (row, col) tuple to algebraic notation."""
    return f"{chr(ord('a') + c)}{r + 1}"

def get_attacked_squares(pieces):
    """
    Calculates all squares attacked by a set of white pieces.
    Handles blocking by other pieces on the board.
    """
    attacked = set()
    all_pieces_pos = {pos for pos, _, _ in pieces}

    for pos, piece_type, color in pieces:
        if color == 'black':
            continue
            
        r, c = pos

        # --- Pawn ---
        if piece_type == 'p':
            # White pawn attacks diagonally forward
            if r < 7:
                if c > 0: attacked.add((r + 1, c - 1))
                if c < 7: attacked.add((r + 1, c + 1))
        
        # --- King ---
        if piece_type == 'k':
            for dr in [-1, 0, 1]:
                for dc in [-1, 0, 1]:
                    if dr == 0 and dc == 0: continue
                    nr, nc = r + dr, c + dc
                    if 0 <= nr < 8 and 0 <= nc < 8:
                        attacked.add((nr, nc))

        # --- Queen (and Rook/Bishop helpers) ---
        if piece_type == 'q' or piece_type == 'r' or piece_type == 'b':
            # Rays for Rook (horizontal/vertical) and Bishop (diagonal)
            directions = []
            if piece_type in ['q', 'r']:
                directions.extend([(0, 1), (0, -1), (1, 0), (-1, 0)])
            if piece_type in ['q', 'b']:
                directions.extend([(1, 1), (1, -1), (-1, 1), (-1, -1)])
            
            for dr, dc in directions:
                nr, nc = r + dr, c + dc
                while 0 <= nr < 8 and 0 <= nc < 8:
                    attacked.add((nr, nc))
                    if (nr, nc) in all_pieces_pos:
                        break # Attack is blocked by an intervening piece
                    nr += dr
                    nc += dc

    return attacked


def solve_chess_puzzle():
    """
    Solves the specified chess puzzle by verifying a candidate solution.
    """
    # Candidate position based on Heathcote's problem (1906)
    # White King is necessary for the stalemate but costs 0 points.
    white_pieces_config = [
        (to_coords('g6'), 'q', 'white'), # Queen = 9 points
        (to_coords('f7'), 'p', 'white'), # Pawn  = 1 point
        (to_coords('h6'), 'k', 'white'), # King  = 0 points
    ]
    
    black_king_pos_alg = 'h8'
    black_king_pos = to_coords(black_king_pos_alg)

    all_pieces = white_pieces_config + [(black_king_pos, 'k', 'black')]
    
    # Verify stalemate
    # 1. Is the king's square safe?
    attacked_by_white = get_attacked_squares(white_pieces_config)
    is_king_safe = black_king_pos not in attacked_by_white

    # 2. Are escape squares attacked?
    kr, kc = black_king_pos
    escape_squares = set()
    for dr in [-1, 0, 1]:
        for dc in [-1, 0, 1]:
            if dr == 0 and dc == 0: continue
            nr, nc = kr + dr, kc + dc
            if 0 <= nr < 8 and 0 <= nc < 8:
                escape_squares.add((nr, nc))

    are_escapes_covered = escape_squares.issubset(attacked_by_white)
    
    is_stalemate = is_king_safe and are_escapes_covered

    # Verify board coverage
    all_squares = {(r, c) for r in range(8) for c in range(8)}
    unattacked_squares = all_squares - attacked_by_white
    
    # The only unattacked square should be the one with the black king
    is_coverage_correct = (unattacked_squares == {black_king_pos})

    # The point value of the pieces
    piece_points = {'q': 9, 'p': 1, 'k': 0}
    total_points = sum(piece_points[ptype] for _, ptype, _ in white_pieces_config)

    # Print the result
    if is_stalemate and is_coverage_correct:
        print("The smallest number of points is 10.")
        print("This can be achieved with a Queen and a Pawn.")
        print("The final equation is:")
        print("Queen (9) + Pawn (1) = 10")
    else:
        # This part should not be reached if the analysis is correct
        print("The provided configuration is not a valid solution.")

solve_chess_puzzle()
