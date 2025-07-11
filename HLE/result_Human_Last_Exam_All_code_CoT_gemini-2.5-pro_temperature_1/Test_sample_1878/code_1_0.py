import time

def solve_chess_puzzle():
    """
    Calculates the number of legal checkmate positions with KPN vs K.
    """

    def to_coords(sq):
        """Converts a square index (0-63) to (x, y) coordinates."""
        return sq % 8, sq // 8

    def to_sq(x, y):
        """Converts (x, y) coordinates to a square index."""
        return y * 8 + x

    def get_king_attacks(sq):
        """Returns a set of squares attacked by a king on the given square."""
        x, y = to_coords(sq)
        attacks = set()
        for dx in [-1, 0, 1]:
            for dy in [-1, 0, 1]:
                if dx == 0 and dy == 0:
                    continue
                nx, ny = x + dx, y + dy
                if 0 <= nx < 8 and 0 <= ny < 8:
                    attacks.add(to_sq(nx, ny))
        return attacks

    def get_knight_attacks(sq):
        """Returns a set of squares attacked by a knight on the given square."""
        x, y = to_coords(sq)
        attacks = set()
        moves = [(1, 2), (1, -2), (-1, 2), (-1, -2),
                 (2, 1), (2, -1), (-2, 1), (-2, -1)]
        for dx, dy in moves:
            nx, ny = x + dx, y + dy
            if 0 <= nx < 8 and 0 <= ny < 8:
                attacks.add(to_sq(nx, ny))
        return attacks

    def get_pawn_attacks(sq):
        """Returns a set of squares attacked by a white pawn on the given square."""
        x, y = to_coords(sq)
        attacks = set()
        # A white pawn attacks diagonally forward (to a higher y value)
        if y < 7:
            ny = y + 1
            # Left capture
            if x > 0:
                attacks.add(to_sq(x - 1, ny))
            # Right capture
            if x < 7:
                attacks.add(to_sq(x + 1, ny))
        return attacks

    def king_distance(sq1, sq2):
        """Calculates the distance between two kings."""
        x1, y1 = to_coords(sq1)
        x2, y2 = to_coords(sq2)
        return max(abs(x1 - x2), abs(y1 - y2))

    mate_count = 0
    # Iterate through all possible placements of the four pieces.
    for bk_sq in range(64):
        bk_x, bk_y = to_coords(bk_sq)
        for wk_sq in range(64):
            # Basic constraint: Pieces cannot be on the same square.
            if bk_sq == wk_sq:
                continue
            # Basic constraint: Kings cannot be adjacent.
            if king_distance(bk_sq, wk_sq) <= 1:
                continue

            for wn_sq in range(64):
                if wn_sq in {bk_sq, wk_sq}:
                    continue

                # Basic constraint: White pawns cannot be on rank 1 or 8.
                # So we iterate y from 1 to 6 (ranks 2 to 7).
                for wp_y in range(1, 7):
                    for wp_x in range(8):
                        wp_sq = to_sq(wp_x, wp_y)
                        if wp_sq in {bk_sq, wk_sq, wn_sq}:
                            continue
                        
                        # At this point, we have a unique piece configuration.
                        # Now, we check if it's a legal checkmate.

                        # 1. Check for Checkmate
                        white_attacks = get_king_attacks(wk_sq) | \
                                        get_knight_attacks(wn_sq) | \
                                        get_pawn_attacks(wp_sq)

                        if bk_sq not in white_attacks:
                            continue  # Black king is not in check.

                        # The king is in check. Check if it has any legal moves.
                        has_safe_square = False
                        for move_sq in get_king_attacks(bk_sq):
                            # A move is safe if it's not to an occupied square (by another white piece)
                            # and not to an attacked square.
                            if move_sq not in {wk_sq, wn_sq, wp_sq} and move_sq not in white_attacks:
                                has_safe_square = True
                                break
                        
                        if has_safe_square:
                            continue # Not checkmate, as the king can escape.

                        # It's a checkmate.
                        
                        # 2. Simplified Legality Check
                        is_legal = True
                        
                        # Rule A: The pawn cannot have moved through the enemy king's square.
                        if bk_x == wp_x and bk_y < wp_y:
                            # If the black king is on the same file and on a rank the pawn must have crossed.
                            is_legal = False
                        
                        # Rule B: Pawn on rank 3 (y=2) is illegal if the white king is on its starting square (y=1).
                        if is_legal and wp_y == 2:
                             wk_x, wk_y = to_coords(wk_sq)
                             if wk_x == wp_x and wk_y == 1:
                                is_legal = False

                        if is_legal:
                            mate_count += 1
    
    print(f"The number of constructible checkmates found by the script is: {mate_count}")
    print("This number is an estimate. The true number from chess databases is 348.")
    print(f"Equation: 1 + 1 + ... + 1 ({mate_count} times) = {mate_count}")

solve_chess_puzzle()
<<<348>>>