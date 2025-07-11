import sys

def solve():
    """
    Calculates the number of legal checkmate positions with White (King, Pawn, Knight)
    against a lone Black King.
    """

    # Helper functions to convert between square index (0-63) and coordinates (row, col)
    def to_coords(sq):
        return sq // 8, sq % 8

    def to_sq(r, c):
        return r * 8 + c

    # Pre-calculate attack sets for each piece from each square to speed up the main loop.
    # A chess board is represented as 0=a1, 1=b1, ..., 63=h8.
    # White pawns move from lower row index to higher row index.

    KING_ATTACKS = {}
    for sq in range(64):
        r, c = to_coords(sq)
        attacks = set()
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0:
                    continue
                nr, nc = r + dr, c + dc
                if 0 <= nr < 8 and 0 <= nc < 8:
                    attacks.add(to_sq(nr, nc))
        KING_ATTACKS[sq] = attacks

    KNIGHT_ATTACKS = {}
    for sq in range(64):
        r, c = to_coords(sq)
        attacks = set()
        moves = [(1, 2), (1, -2), (-1, 2), (-1, -2), (2, 1), (2, -1), (-2, 1), (-2, -1)]
        for dr, dc in moves:
            nr, nc = r + dr, c + dc
            if 0 <= nr < 8 and 0 <= nc < 8:
                attacks.add(to_sq(nr, nc))
        KNIGHT_ATTACKS[sq] = attacks

    # White pawns attack diagonally forward (to a higher row index).
    PAWN_ATTACKS = {}
    for sq in range(64):
        r, c = to_coords(sq)
        attacks = set()
        nr = r + 1  # White pawns move "up" the board from row 0 to 7
        if 0 <= nr < 8:
            if c > 0:
                attacks.add(to_sq(nr, c - 1))
            if c < 7:
                attacks.add(to_sq(nr, c + 1))
        PAWN_ATTACKS[sq] = attacks

    mate_count = 0
    # Iterate through all possible placements of the four pieces.
    for bk_sq in range(64):
        for wk_sq in range(64):
            if bk_sq == wk_sq:
                continue
            
            # Legality Check 1: Kings cannot be adjacent. This also ensures the
            # side not to move (White) is not in check.
            if wk_sq in KING_ATTACKS[bk_sq]:
                continue

            for wn_sq in range(64):
                if wn_sq == bk_sq or wn_sq == wk_sq:
                    continue

                # Legality Check 2: Pawns cannot be on the 1st or 8th rank.
                # We iterate the pawn only on valid ranks (2-7, which are rows 1-6).
                for wp_sq_row in range(1, 7):
                    for wp_sq_col in range(8):
                        wp_sq = to_sq(wp_sq_row, wp_sq_col)

                        if wp_sq == bk_sq or wp_sq == wk_sq or wp_sq == wn_sq:
                            continue

                        # At this point, we have a unique and legal placement.
                        # Now, check for checkmate.

                        # 1. Calculate all squares attacked by White's pieces.
                        white_attacks = (KING_ATTACKS[wk_sq] |
                                         KNIGHT_ATTACKS[wn_sq] |
                                         PAWN_ATTACKS[wp_sq])

                        # 2. Check if the Black King is in check. If not, it's not mate.
                        if bk_sq not in white_attacks:
                            continue

                        # 3. Check if the Black King has any escape squares.
                        # An escape exists if the king can move to an adjacent square
                        # that is NOT attacked by White.
                        is_mate = True
                        for move_sq in KING_ATTACKS[bk_sq]:
                            if move_sq not in white_attacks:
                                is_mate = False
                                break
                        
                        if is_mate:
                            mate_count += 1
    
    # The final result of the calculation.
    print(mate_count)

solve()
<<<356>>>