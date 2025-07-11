import collections

def solve_chess_checkmates():
    """
    Calculates the number of legal checkmate positions where White has a King,
    Pawn, and Knight, and Black has only a King.
    """

    # Pre-calculate attack sets for each piece from each square for efficiency.
    # This avoids recalculating these moves millions of times in the main loop.
    KING_ATTACKS = collections.defaultdict(set)
    KNIGHT_ATTACKS = collections.defaultdict(set)
    PAWN_ATTACKS = collections.defaultdict(set)  # For a white pawn

    for r_start in range(8):
        for f_start in range(8):
            pos_start = r_start * 8 + f_start

            # King moves to all adjacent squares
            for dr in [-1, 0, 1]:
                for df in [-1, 0, 1]:
                    if dr == 0 and df == 0:
                        continue
                    r_end, f_end = r_start + dr, f_start + df
                    if 0 <= r_end < 8 and 0 <= f_end < 8:
                        KING_ATTACKS[pos_start].add(r_end * 8 + f_end)

            # Knight 'L' shaped moves
            for dr, df in [(1, 2), (1, -2), (-1, 2), (-1, -2),
                           (2, 1), (2, -1), (-2, 1), (-2, -1)]:
                r_end, f_end = r_start + dr, f_start + df
                if 0 <= r_end < 8 and 0 <= f_end < 8:
                    KNIGHT_ATTACKS[pos_start].add(r_end * 8 + f_end)

            # White Pawn diagonal attacks
            if r_start < 7:  # Pawns on the 8th rank (rank 7) cannot attack
                for df in [-1, 1]:
                    r_end, f_end = r_start + 1, f_start + df
                    if 0 <= f_end < 8:
                        PAWN_ATTACKS[pos_start].add(r_end * 8 + f_end)

    def is_legal_checkmate(wk_pos, bk_pos, wn_pos, wp_pos):
        """
        Checks if a given position is a legal checkmate.
        """
        # A position is illegal if the kings are adjacent. This also covers
        # the rule that the side not to move (White) cannot be in check.
        if wk_pos in KING_ATTACKS[bk_pos]:
            return False

        # Determine all squares attacked by White's pieces.
        white_attacks = (KING_ATTACKS[wk_pos] |
                         KNIGHT_ATTACKS[wn_pos] |
                         PAWN_ATTACKS[wp_pos])

        # Condition 1: The Black King must be in check.
        if bk_pos not in white_attacks:
            return False

        # Condition 2: The Black King must have no legal moves.
        # Get all potential escape squares for the Black King.
        bk_escape_squares = KING_ATTACKS[bk_pos]
        
        # The King cannot move to a square occupied by a white piece.
        white_piece_positions = {wk_pos, wn_pos, wp_pos}
        bk_escape_squares -= white_piece_positions

        # Every remaining escape square must be attacked by White.
        # If we find even one safe square, it's not checkmate.
        for move in bk_escape_squares:
            if move not in white_attacks:
                return False
        
        # If the King is in check and has no safe escape squares, it's checkmate.
        return True

    checkmate_count = 0
    # Iterate through all possible, unique piece placements.
    for bk_pos in range(64):
        for wk_pos in range(64):
            if wk_pos == bk_pos: continue
            for wn_pos in range(64):
                if wn_pos == bk_pos or wn_pos == wk_pos: continue
                # Pawns cannot be on the 1st (rank 0) or 8th (rank 7) ranks.
                for wp_pos in range(8, 56):
                    if wp_pos == bk_pos or wp_pos == wk_pos or wp_pos == wn_pos: continue
                    
                    if is_legal_checkmate(wk_pos, bk_pos, wn_pos, wp_pos):
                        checkmate_count += 1
    
    # The final result is the total count of such positions.
    print(f"The number of constructible checkmate positions is: {checkmate_count}")

solve_chess_checkmates()