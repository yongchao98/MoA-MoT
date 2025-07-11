import itertools

def solve_chess_checkmates():
    """
    Calculates the number of legal checkmate positions where White has a King,
    Knight, and Pawn, and Black has only a King.
    """

    # Helper function to check if a square is attacked by a knight
    def is_attacked_by_knight(target_sq, knight_sq):
        target_r, target_c = target_sq // 8, target_sq % 8
        knight_r, knight_c = knight_sq // 8, knight_sq % 8
        dr = abs(target_r - knight_r)
        dc = abs(target_c - knight_c)
        return (dr == 1 and dc == 2) or (dr == 2 and dc == 1)

    # Helper function to check if a square is attacked by a white pawn
    def is_attacked_by_pawn(target_sq, pawn_sq):
        pawn_r, pawn_c = pawn_sq // 8, pawn_sq % 8
        target_r, target_c = target_sq // 8, target_sq % 8
        # A white pawn attacks one rank forward (r+1) and one file to the side
        return target_r == pawn_r + 1 and abs(target_c - pawn_c) == 1

    # Helper function to get a set of all squares attacked by white pieces
    def get_all_attacked_squares(wk, wn, wp):
        attacked = set()
        
        # King attacks
        wk_r, wk_c = wk // 8, wk % 8
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0: continue
                r, c = wk_r + dr, wk_c + dc
                if 0 <= r < 8 and 0 <= c < 8:
                    attacked.add(r * 8 + c)
        
        # Knight attacks
        wn_r, wn_c = wn // 8, wn % 8
        for dr, dc in [(1, 2), (1, -2), (-1, 2), (-1, -2), (2, 1), (2, -1), (-2, 1), (-2, -1)]:
            r, c = wn_r + dr, wn_c + dc
            if 0 <= r < 8 and 0 <= c < 8:
                attacked.add(r * 8 + c)
                
        # Pawn attacks
        wp_r, wp_c = wp // 8, wp % 8
        if wp_r < 7: # A pawn on the 8th rank is promoted, not a pawn
            if wp_c > 0: attacked.add((wp_r + 1) * 8 + (wp_c - 1))
            if wp_c < 7: attacked.add((wp_r + 1) * 8 + (wp_c + 1))
            
        return attacked

    knight_only_mates = 0
    pawn_only_mates = 0
    double_check_mates = 0

    # Iterate through all permutations of 4 pieces on 64 squares.
    # wk=white king, wn=white knight, wp=white pawn, bk=black king
    # This is a long-running process.
    for wk, wn, wp, bk in itertools.permutations(range(64), 4):
        
        # Rule 1: Pawn cannot be on the 1st or 8th rank.
        if wp // 8 in {0, 7}:
            continue
            
        # Rule 2: Kings cannot be on adjacent squares.
        wk_r, wk_c = wk // 8, wk % 8
        bk_r, bk_c = bk // 8, bk % 8
        if max(abs(wk_r - bk_r), abs(wk_c - bk_c)) <= 1:
            continue

        # --- Checkmate Evaluation ---
        
        # Pre-calculate all squares attacked by White for efficiency
        white_attacks = get_all_attacked_squares(wk, wn, wp)
        
        # Condition 1: The Black King must be in check.
        if bk not in white_attacks:
            continue
            
        # Condition 2: The Black King must have no legal moves.
        has_escape = False
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0: continue
                escape_r, escape_c = bk_r + dr, bk_c + dc
                
                if 0 <= escape_r < 8 and 0 <= escape_c < 8:
                    escape_sq = escape_r * 8 + escape_c
                    if escape_sq not in white_attacks:
                        has_escape = True
                        break
            if has_escape:
                break
        
        if has_escape:
            continue
            
        # If we reach here, it's a valid checkmate. Now categorize it.
        is_knight_check = is_attacked_by_knight(bk, wn)
        is_pawn_check = is_attacked_by_pawn(bk, wp)

        if is_knight_check and is_pawn_check:
            double_check_mates += 1
        elif is_knight_check:
            knight_only_mates += 1
        elif is_pawn_check:
            pawn_only_mates += 1
            
    total_mates = knight_only_mates + pawn_only_mates + double_check_mates
    
    print(f"Found {knight_only_mates} checkmates by Knight + {pawn_only_mates} checkmates by Pawn + {double_check_mates} checkmates by double check = {total_mates} total checkmates")

if __name__ == '__main__':
    solve_chess_checkmates()