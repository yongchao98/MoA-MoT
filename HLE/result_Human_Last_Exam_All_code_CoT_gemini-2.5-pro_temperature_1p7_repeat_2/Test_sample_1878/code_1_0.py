import itertools

def is_valid_sq(r, c):
    """Checks if a square (row, col) is on the 8x8 board."""
    return 0 <= r < 8 and 0 <= c < 8

def get_white_attacks(wk_pos, wn_pos, wp_pos):
    """
    Calculates the set of all squares attacked by White's pieces.
    - wk_pos: (row, col) of White King
    - wn_pos: (row, col) of White Knight
    - wp_pos: (row, col) of White Pawn
    
    Board Orientation: row 0 is Black's back rank (rank 8), row 7 is White's (rank 1).
    White pawns move from higher row index to lower row index.
    """
    attacked_squares = set()
    wk_r, wk_c = wk_pos
    wn_r, wn_c = wn_pos
    wp_r, wp_c = wp_pos

    # 1. King attacks
    for dr in [-1, 0, 1]:
        for dc in [-1, 0, 1]:
            if dr == 0 and dc == 0:
                continue
            nr, nc = wk_r + dr, wk_c + dc
            if is_valid_sq(nr, nc):
                attacked_squares.add((nr, nc))

    # 2. Knight attacks
    knight_moves = [
        (1, 2), (1, -2), (-1, 2), (-1, -2),
        (2, 1), (2, -1), (-2, 1), (-2, -1)
    ]
    for dr, dc in knight_moves:
        nr, nc = wn_r + dr, wn_c + dc
        if is_valid_sq(nr, nc):
            attacked_squares.add((nr, nc))

    # 3. Pawn attacks (diagonal forward)
    for dc in [-1, 1]:
        nr, nc = wp_r - 1, wp_c + dc
        if is_valid_sq(nr, nc):
            attacked_squares.add((nr, nc))
            
    return attacked_squares

def find_checkmates():
    """
    Finds the number of legal checkmate positions with WK, WN, WP vs BK.
    A position is 'legal' if:
    1. Kings are not adjacent.
    2. White's pawn is not on the 1st or 8th rank.
    """
    mate_count = 0
    all_squares = range(64)

    # This loop iterates through 15,249,024 permutations. It may take a few minutes to run.
    for wk_sq, wn_sq, wp_sq, bk_sq in itertools.permutations(all_squares, 4):
        wk_pos = divmod(wk_sq, 8)
        wn_pos = divmod(wn_sq, 8)
        wp_pos = divmod(wp_sq, 8)
        bk_pos = divmod(bk_sq, 8)

        # ========== 1. Legality Check ==========

        # White Pawn cannot be on the 1st (row 7) or 8th (row 0) rank.
        if wp_pos[0] == 0 or wp_pos[0] == 7:
            continue

        # Kings cannot be on adjacent squares.
        # This also ensures White's king is not in check.
        if max(abs(wk_pos[0] - bk_pos[0]), abs(wk_pos[1] - bk_pos[1])) <= 1:
            continue

        # ========== 2. Checkmate Validation ==========

        white_attacked_squares = get_white_attacks(wk_pos, wn_pos, wp_pos)

        # The Black King must be in check.
        if bk_pos not in white_attacked_squares:
            continue

        # The Black King must have no legal moves.
        has_legal_move = False
        bk_r, bk_c = bk_pos
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0:
                    continue
                
                escape_r, escape_c = bk_r + dr, bk_c + dc
                
                # An escape is possible if the destination square is on the board
                # AND it is not under attack by any white piece.
                if is_valid_sq(escape_r, escape_c) and (escape_r, escape_c) not in white_attacked_squares:
                    has_legal_move = True
                    break
            if has_legal_move:
                break
        
        # If the king is in check and has no legal moves, it's checkmate.
        if not has_legal_move:
            mate_count += 1
            
    print(f"The total number of legal checkmate positions is: {mate_count}")

# Calling the function to start the calculation.
find_checkmates()