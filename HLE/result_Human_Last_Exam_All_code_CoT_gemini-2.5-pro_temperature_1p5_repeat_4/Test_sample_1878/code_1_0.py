def solve_chess_checkmates():
    """
    Calculates the number of legal checkmate positions with White (King, Pawn, Knight)
    vs. Black (King).
    """
    mate_count = 0
    squares = range(64)

    # --- Pre-computation for speed ---
    # For each square, pre-calculate its (row, col) coordinates
    coords = {sq: (sq // 8, sq % 8) for sq in squares}

    # For each square, pre-calculate where a knight could move from it
    knight_moves = {}
    for r_start in range(8):
        for c_start in range(8):
            sq = r_start * 8 + c_start
            moves = set()
            for dr, dc in [(1, 2), (1, -2), (-1, 2), (-1, -2),
                           (2, 1), (2, -1), (-2, 1), (-2, -1)]:
                r_end, c_end = r_start + dr, c_start + dc
                if 0 <= r_end < 8 and 0 <= c_end < 8:
                    moves.add(r_end * 8 + c_end)
            knight_moves[sq] = moves

    # For each square, pre-calculate where a king could move from it
    king_moves = {}
    for r_start in range(8):
        for c_start in range(8):
            sq = r_start * 8 + c_start
            moves = set()
            for dr in [-1, 0, 1]:
                for dc in [-1, 0, 1]:
                    if dr == 0 and dc == 0:
                        continue
                    r_end, c_end = r_start + dr, c_start + dc
                    if 0 <= r_end < 8 and 0 <= c_end < 8:
                        moves.add(r_end * 8 + c_end)
            king_moves[sq] = moves

    # For each square, pre-calculate where a white pawn would attack from it
    pawn_attacks = {}
    for r_start in range(8):
        for c_start in range(8):
            sq = r_start * 8 + c_start
            attacks = set()
            r_end = r_start + 1  # White pawns move "up" the board
            if 0 <= r_end < 8:
                # Attack diagonally
                if c_start > 0:
                    attacks.add(r_end * 8 + (c_start - 1))
                if c_start < 7:
                    attacks.add(r_end * 8 + (c_start + 1))
            pawn_attacks[sq] = attacks

    # --- Main Loop: Iterate through all unique piece placements ---
    for bk_sq in squares:
        bk_r, bk_c = coords[bk_sq]
        for wk_sq in squares:
            if wk_sq == bk_sq: continue
            
            # Legality Check 1: Kings must not be adjacent.
            wk_r, wk_c = coords[wk_sq]
            if max(abs(bk_r - wk_r), abs(bk_c - wk_c)) <= 1:
                continue

            for wn_sq in squares:
                if wn_sq == bk_sq or wn_sq == wk_sq: continue

                for wp_sq in squares:
                    if wp_sq in (bk_sq, wk_sq, wn_sq): continue

                    # Legality Check 2: White Pawn cannot be on the 1st or 8th rank.
                    if coords[wp_sq][0] in (0, 7):
                        continue

                    # --- Checkmate Verification ---
                    # 1. Determine all squares attacked by White.
                    attacked_by_white = king_moves[wk_sq] | knight_moves[wn_sq] | pawn_attacks[wp_sq]
                    
                    # 2. The Black King must be in check.
                    if bk_sq not in attacked_by_white:
                        continue
                        
                    # 3. The Black King must have no safe square to move to.
                    has_escape_square = False
                    for move_sq in king_moves[bk_sq]:
                        if move_sq not in attacked_by_white:
                            has_escape_square = True
                            break
                    
                    # 4. If in check with no escape, it's checkmate.
                    if not has_escape_square:
                        mate_count += 1
                        
    print(f"Total legal checkmate positions found: {mate_count}")

solve_chess_checkmates()