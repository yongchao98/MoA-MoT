import itertools

def get_row_col(pos):
    """Converts a 0-63 square index to a (row, col) tuple."""
    return pos // 8, pos % 8

def get_pos(row, col):
    """Converts a (row, col) tuple to a 0-63 square index."""
    return row * 8 + col

def get_king_attacks(pos):
    """Returns a set of squares attacked by a king at a given position."""
    attacks = set()
    r, c = get_row_col(pos)
    for dr in [-1, 0, 1]:
        for dc in [-1, 0, 1]:
            if dr == 0 and dc == 0:
                continue
            nr, nc = r + dr, c + dc
            if 0 <= nr < 8 and 0 <= nc < 8:
                attacks.add(get_pos(nr, nc))
    return attacks

def get_knight_attacks(pos):
    """Returns a set of squares attacked by a knight at a given position."""
    attacks = set()
    r, c = get_row_col(pos)
    moves = [(-2, -1), (-2, 1), (-1, -2), (-1, 2),
             (1, -2), (1, 2), (2, -1), (2, 1)]
    for dr, dc in moves:
        nr, nc = r + dr, c + dc
        if 0 <= nr < 8 and 0 <= nc < 8:
            attacks.add(get_pos(nr, nc))
    return attacks

def get_pawn_attacks(pos):
    """Returns a set of squares attacked by a white pawn at a given position."""
    attacks = set()
    r, c = get_row_col(pos)
    # White pawn moves from rank 2 (row 1) to 7 (row 6) and attacks forward-diagonally.
    if r < 7:
        nr = r + 1
        for dc in [-1, 1]:
            nc = c + dc
            if 0 <= nc < 8:
                attacks.add(get_pos(nr, nc))
    return attacks

def is_kings_adjacent(k1_pos, k2_pos):
    """Checks if two kings are on adjacent squares."""
    r1, c1 = get_row_col(k1_pos)
    r2, c2 = get_row_col(k2_pos)
    return abs(r1 - r2) <= 1 and abs(c1 - c2) <= 1

def solve_chess_puzzle():
    """
    Calculates the number of legal checkmate positions where White has K, N, P
    and Black has only K.
    """
    checkmate_count = 0
    squares = range(64)

    # Pre-computation for performance
    KING_ATTACKS = [get_king_attacks(i) for i in squares]
    KNIGHT_ATTACKS = [get_knight_attacks(i) for i in squares]
    PAWN_ATTACKS = [get_pawn_attacks(i) for i in squares]
    # Pawns cannot legally exist on the 1st or 8th rank.
    PAWN_INVALID_SQUARES = {i for i in squares if get_row_col(i)[0] in [0, 7]}

    # Iterate through all unique placements of the 4 pieces
    # (White King, Black King, White Knight, White Pawn)
    for wk_pos, bk_pos, wn_pos, wp_pos in itertools.permutations(squares, 4):

        # --- Filter 1: Pawn on illegal rank ---
        if wp_pos in PAWN_INVALID_SQUARES:
            continue

        # --- Filter 2: Kings are adjacent ---
        if is_kings_adjacent(wk_pos, bk_pos):
            continue

        # --- Check for Checkmate ---

        # The attacking pieces' moves (pawn or knight must deliver check)
        checking_attacks = KNIGHT_ATTACKS[wn_pos].union(PAWN_ATTACKS[wp_pos])

        # Condition 1: Is the Black King in check?
        if bk_pos not in checking_attacks:
            continue

        # Condition 2: Can the Black King escape?
        # An escape square must not be attacked or occupied by White.

        white_pieces_pos = {wk_pos, wn_pos, wp_pos}
        # All squares controlled by White
        all_white_controlled_squares = checking_attacks.union(KING_ATTACKS[wk_pos])
        # All squares the Black King cannot move to
        forbidden_escape_squares = all_white_controlled_squares.union(white_pieces_pos)

        can_escape = False
        # Check all adjacent squares for the Black King
        for escape_pos in KING_ATTACKS[bk_pos]:
            if escape_pos not in forbidden_escape_squares:
                can_escape = True
                break
        
        # If the king is in check and cannot escape, it's checkmate.
        if not can_escape:
            checkmate_count += 1
            
    print("Final result:")
    print("The total number of legal checkmate positions is:")
    print(checkmate_count)

# Execute the main function
solve_chess_puzzle()