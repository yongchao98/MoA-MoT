import itertools

def count_checkmate_positions():
    """
    Calculates the number of legal checkmate positions with White (K, P, N) vs Black (K).

    A position is considered a legal checkmate if:
    1. All pieces are on distinct squares.
    2. The white pawn is on ranks 2-7.
    3. The kings are not on adjacent squares.
    4. The black king is in check by the knight or pawn.
    5. The black king has no legal moves (all adjacent squares are attacked by white).
    """

    # Pre-calculate knight moves for all 64 squares for O(1) lookup later.
    # knight_attacks[square_index] = {set of attacked squares}
    knight_attacks = [set() for _ in range(64)]
    for r1 in range(8):
        for c1 in range(8):
            s1 = r1 * 8 + c1
            for dr, dc in [(1, 2), (1, -2), (-1, 2), (-1, -2),
                           (2, 1), (2, -1), (-2, 1), (-2, -1)]:
                r2, c2 = r1 + dr, c1 + dc
                if 0 <= r2 < 8 and 0 <= c2 < 8:
                    s2 = r2 * 8 + c2
                    knight_attacks[s1].add(s2)

    # Pre-calculate king moves (attacked squares).
    king_attacks = [set() for _ in range(64)]
    for r1 in range(8):
        for c1 in range(8):
            s1 = r1 * 8 + c1
            for dr in [-1, 0, 1]:
                for dc in [-1, 0, 1]:
                    if dr == 0 and dc == 0:
                        continue
                    r2, c2 = r1 + dr, c1 + dc
                    if 0 <= r2 < 8 and 0 <= c2 < 8:
                        s2 = r2 * 8 + c2
                        king_attacks[s1].add(s2)

    # Pre-calculate pawn attacks. White pawns move from row r to r+1.
    # Pawns can only deliver check from rows 1-6 (ranks 2-7).
    pawn_attacks = [set() for _ in range(64)]
    for r1 in range(1, 7): 
        for c1 in range(8):
            s1 = r1 * 8 + c1
            if c1 > 0: # Left capture
                pawn_attacks[s1].add((r1 + 1) * 8 + (c1 - 1))
            if c1 < 7: # Right capture
                pawn_attacks[s1].add((r1 + 1) * 8 + (c1 + 1))

    mate_count = 0
    
    # Iterate through all permutations of 4 pieces on 64 squares.
    # p = (white_king_sq, white_knight_sq, white_pawn_sq, black_king_sq)
    all_squares = range(64)
    for p in itertools.permutations(all_squares, 4):
        wk_sq, wn_sq, wp_sq, bk_sq = p

        # Filter 1: White pawn must be on a legal rank (2-7, i.e., row 1-6).
        wp_r = wp_sq // 8
        if not (1 <= wp_r <= 6):
            continue

        # Filter 2: Kings cannot be adjacent.
        if bk_sq in king_attacks[wk_sq]:
            continue

        # Filter 3: Black king must be in check by the knight or pawn.
        is_in_check = (bk_sq in knight_attacks[wn_sq]) or \
                      (bk_sq in pawn_attacks[wp_sq])
        
        if not is_in_check:
            continue
            
        # Filter 4: Black king must have no legal moves (checkmate).
        # Union of all squares attacked by white pieces.
        white_controlled_squares = (king_attacks[wk_sq] | 
                                    knight_attacks[wn_sq] | 
                                    pawn_attacks[wp_sq])
        
        has_escape_square = False
        # Check all squares the black king could move to.
        for dest_sq in king_attacks[bk_sq]:
            if dest_sq not in white_controlled_squares:
                has_escape_square = True
                break
        
        if has_escape_square:
            continue

        # If all checks pass, this is a valid checkmate position.
        mate_count += 1
        
    print(f"Total number of checkmate positions: {mate_count}")

# Note: This script is computationally intensive and may take a few minutes to run.
if __name__ == '__main__':
    count_checkmate_positions()