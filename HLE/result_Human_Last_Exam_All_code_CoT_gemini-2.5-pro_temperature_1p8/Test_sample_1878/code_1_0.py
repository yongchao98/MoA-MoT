import itertools

def to_coords(square_index):
    """Converts a 0-63 square index to (rank, file) coordinates."""
    return (square_index // 8, square_index % 8)

def get_king_attacks(r, f):
    """Returns a set of squares attacked by a king at (rank, file)."""
    attacks = set()
    for dr in [-1, 0, 1]:
        for df in [-1, 0, 1]:
            if dr == 0 and df == 0:
                continue
            nr, nf = r + dr, f + df
            if 0 <= nr <= 7 and 0 <= nf <= 7:
                attacks.add((nr, nf))
    return attacks

def get_knight_attacks(r, f):
    """Returns a set of squares attacked by a knight at (rank, file)."""
    attacks = set()
    moves = [(2, 1), (2, -1), (-2, 1), (-2, -1),
             (1, 2), (1, -2), (-1, 2), (-1, -2)]
    for dr, df in moves:
        nr, nf = r + dr, f + df
        if 0 <= nr <= 7 and 0 <= nf <= 7:
            attacks.add((nr, nf))
    return attacks

def get_pawn_attacks(r, f):
    """Returns a set of squares attacked by a white pawn at (rank, file)."""
    attacks = set()
    # A white pawn attacks diagonally forward (increasing rank).
    for df in [-1, 1]:
        nr, nf = r + 1, f + df
        if 0 <= nr <= 7 and 0 <= nf <= 7:
            attacks.add((nr, nf))
    return attacks

def solve_chess_puzzle():
    """
    Calculates the number of legal checkmate positions with WK, WN, WP vs BK.
    """
    squares = range(64)
    mate_positions_count = 0

    # Generate all unique placements of the 4 pieces on 64 squares.
    for piece_squares in itertools.permutations(squares, 4):
        wk_sq, wn_sq, wp_sq, bk_sq = piece_squares

        wk_pos = to_coords(wk_sq)
        wn_pos = to_coords(wn_sq)
        wp_pos = to_coords(wp_sq)
        bk_pos = to_coords(bk_sq)

        # Filter 1: White Pawn cannot be on rank 1 (index 0) or 8 (index 7).
        if wp_pos[0] in [0, 7]:
            continue

        # Filter 2: Kings cannot be on adjacent squares.
        # This checks if the maximum distance in rank or file is 1 or less.
        if max(abs(wk_pos[0] - bk_pos[0]), abs(wk_pos[1] - bk_pos[1])) <= 1:
            continue

        # --- Checkmate Verification ---
        
        # 1. Calculate all squares attacked by White's pieces.
        white_king_attacks = get_king_attacks(wk_pos[0], wk_pos[1])
        white_knight_attacks = get_knight_attacks(wn_pos[0], wn_pos[1])
        white_pawn_attacks = get_pawn_attacks(wp_pos[0], wp_pos[1])
        all_white_attacks = white_king_attacks.union(white_knight_attacks).union(white_pawn_attacks)

        # 2. Check if the Black King is currently in check.
        if bk_pos not in all_white_attacks:
            continue

        # 3. Check if the Black King has any legal moves (escape squares).
        has_escape_square = False
        white_piece_positions = {wk_pos, wn_pos, wp_pos}
        
        for move in get_king_attacks(bk_pos[0], bk_pos[1]):
            # An escape square must not be occupied by a white piece AND not be attacked.
            if move not in white_piece_positions and move not in all_white_attacks:
                has_escape_square = True
                break  # Found an escape, so it's not mate.

        # If the king is in check and has no escape squares, it's checkmate.
        if not has_escape_square:
            mate_positions_count += 1
    
    print(mate_positions_count)

# Run the calculation.
solve_chess_puzzle()