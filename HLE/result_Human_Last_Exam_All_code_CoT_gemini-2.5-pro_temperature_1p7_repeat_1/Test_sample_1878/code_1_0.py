import itertools

def generate_attack_sets():
    """
    Pre-computes and returns attack sets for kings, knights, and white pawns
    for every square on the board (0-63).
    """
    king_attacks = {}
    knight_attacks = {}
    pawn_attacks = {}  # For White Pawns

    for sq in range(64):
        r, c = sq // 8, sq % 8

        # King attacks
        k_attacks = set()
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0:
                    continue
                nr, nc = r + dr, c + dc
                if 0 <= nr < 8 and 0 <= nc < 8:
                    k_attacks.add(nr * 8 + nc)
        king_attacks[sq] = k_attacks

        # Knight attacks
        n_attacks = set()
        for dr, dc in [(1, 2), (1, -2), (-1, 2), (-1, -2), (2, 1), (2, -1), (-2, 1), (-2, -1)]:
            nr, nc = r + dr, c + dc
            if 0 <= nr < 8 and 0 <= nc < 8:
                n_attacks.add(nr * 8 + nc)
        knight_attacks[sq] = n_attacks

        # White Pawn attacks (moves "up" the board, from rank 1 to 8)
        p_attacks = set()
        nr = r + 1  # Pawn on rank `r` attacks rank `r+1`
        if nr < 8:
            if c > 0: p_attacks.add(nr * 8 + c - 1)
            if c < 7: p_attacks.add(nr * 8 + c + 1)
        pawn_attacks[sq] = p_attacks
        
    return king_attacks, knight_attacks, pawn_attacks

def find_checkmates():
    """
    Iterates through all legal piece placements to find and categorize checkmates.
    """
    KING_ATTACKS, KNIGHT_ATTACKS, PAWN_ATTACKS = generate_attack_sets()
    
    counts = {
        "knight_mate": 0,
        "pawn_mate": 0,
        "double_check_mate": 0,
    }

    # Iterate through all unique placements of the 4 pieces on 64 squares.
    # p = (White King, Black King, White Knight, White Pawn)
    for p in itertools.permutations(range(64), 4):
        wk_pos, bk_pos, wn_pos, wp_pos = p

        # --- 1. Legality Checks ---
        
        # White Pawn cannot be on the 1st (rank 0) or 8th (rank 7) rank.
        wp_row = wp_pos // 8
        if wp_row == 0 or wp_row == 7:
            continue
            
        # Kings cannot be on adjacent squares.
        if bk_pos in KING_ATTACKS[wk_pos]:
            continue

        # --- 2. Checkmate Verification ---
        
        # Get attack sets for the current White piece positions
        wk_controlled = KING_ATTACKS[wk_pos]
        wn_controlled = KNIGHT_ATTACKS[wn_pos]
        wp_controlled = PAWN_ATTACKS[wp_pos]
        
        # a) Is the Black King in check?
        checking_pieces = []
        if bk_pos in wn_controlled:
            checking_pieces.append(wn_pos)
        if bk_pos in wp_controlled:
            checking_pieces.append(wp_pos)
            
        if not checking_pieces:
            continue  # Not in check, so not a checkmate.

        # b) Does the Black King have any legal escape squares?
        all_white_controlled_squares = wk_controlled | wn_controlled | wp_controlled
        
        has_escape = False
        for move_sq in KING_ATTACKS[bk_pos]:
            if move_sq not in all_white_controlled_squares:
                has_escape = True
                break
        
        if has_escape:
            continue  # Not checkmate, the king can move away.

        # c) Can the check be legally removed by capture?
        #    If it's a double check, capture is impossible, so it's mate.
        if len(checking_pieces) > 1:
            counts["double_check_mate"] += 1
            continue

        # If it's a single check, see if the Black King can capture the checker.
        checker_pos = checking_pieces[0]
        
        # The king can capture if the checker's square is not controlled by other white pieces.
        # "Other" pieces' attacks:
        if checker_pos == wn_pos:  # The Knight is the checker
            attacks_on_checker_sq = wk_controlled | wp_controlled
        else:  # The Pawn is the checker
            attacks_on_checker_sq = wk_controlled | wn_controlled
        
        if checker_pos not in attacks_on_checker_sq:
            # The checker can be captured, so it's not a mate.
            continue
        
        # If we reach here, it's a confirmed checkmate.
        # Let's categorize and count it.
        if checker_pos == wn_pos:
            counts["knight_mate"] += 1
        else: # checker_pos == wp_pos
            counts["pawn_mate"] += 1
    
    # --- 3. Output the Final Result ---
    total_mates = sum(counts.values())
    
    print(f"Number of checkmates delivered by the Knight: {counts['knight_mate']}")
    print(f"Number of checkmates delivered by the Pawn: {counts['pawn_mate']}")
    print(f"Number of checkmates from a Double Check: {counts['double_check_mate']}")
    print("---")
    print("Final Equation:")
    print(f"{counts['knight_mate']} (Knight Mates) + {counts['pawn_mate']} (Pawn Mates) + {counts['double_check_mate']} (Double Check Mates) = {total_mates} (Total)")

    return total_mates

if __name__ == "__main__":
    total = find_checkmates()
    print(f"\nTotal unique checkmate positions found: {total}")
    print(f"<<<{total}>>>")
