import sys

def solve_chess_puzzle():
    """
    This script calculates the number of legal checkmate positions where
    White (King, Pawn, Knight) mates Black (King).
    """

    # Helper function to convert a square index (0-63) to (rank, file)
    def to_coords(sq):
        # rank 0-7, file 0-7
        return sq // 8, sq % 8

    # Helper function to check if two squares are adjacent for kings
    def are_kings_adjacent(sq1, sq2):
        r1, f1 = to_coords(sq1)
        r2, f2 = to_coords(sq2)
        return max(abs(r1 - r2), abs(f1 - f2)) == 1

    # Helper function to check if a square is attacked by any white piece
    def is_square_attacked_by_white(target_sq, wk_sq, wn_sq, wp_sq):
        tr, tf = to_coords(target_sq)

        # 1. Attacked by White King
        if are_kings_adjacent(target_sq, wk_sq):
            return True

        # 2. Attacked by White Knight
        wnr, wnf = to_coords(wn_sq)
        dr, df = abs(tr - wnr), abs(tf - wnf)
        if (dr == 1 and df == 2) or (dr == 2 and df == 1):
            return True

        # 3. Attacked by White Pawn (moving up from rank 0 to 7)
        wpr, wpf = to_coords(wp_sq)
        if tr == wpr + 1 and abs(tf - wpf) == 1:
            return True

        return False

    # Initialize counters for different types of checkmates
    knight_mates = 0
    pawn_mates = 0
    double_check_mates = 0

    # Main loops to iterate through all piece placements
    for bk_sq in range(64):  # Black King
        for wk_sq in range(64):  # White King
            if wk_sq == bk_sq: continue
            if are_kings_adjacent(wk_sq, bk_sq): continue

            for wp_sq in range(8, 56):  # White Pawn (ranks 2-7)
                if wp_sq in {bk_sq, wk_sq}: continue

                # Retrograde legality check: Pawn cannot be on the same file
                # in front of the enemy king, as this position is unreachable.
                wpr, wpf = to_coords(wp_sq)
                bkr, bkf = to_coords(bk_sq)
                if wpf == bkf and bkr == wpr + 1:
                    continue

                for wn_sq in range(64):  # White Knight
                    if wn_sq in {bk_sq, wk_sq, wp_sq}: continue

                    # At this point, we have a unique and legal placement of pieces.
                    # Now, check for checkmate.

                    # Condition 1: Is Black King in check?
                    if not is_square_attacked_by_white(bk_sq, wk_sq, wn_sq, wp_sq):
                        continue

                    # Condition 2: Does Black King have any legal moves?
                    has_safe_escape = False
                    for dr in [-1, 0, 1]:
                        for df in [-1, 0, 1]:
                            if dr == 0 and df == 0: continue
                            
                            nr, nf = bkr + dr, bkf + df

                            if 0 <= nr < 8 and 0 <= nf < 8:
                                escape_sq = nr * 8 + nf
                                # Check if escape square is safe from attack
                                if not is_square_attacked_by_white(escape_sq, wk_sq, wn_sq, wp_sq):
                                    has_safe_escape = True
                                    break
                        if has_safe_escape: break
                    
                    # If king is in check and has no safe escape, it's checkmate.
                    if not has_safe_escape:
                        # Determine which piece(s) delivered the check for categorization.
                        check_by_knight = False
                        check_by_pawn = False
                        
                        # Knight check?
                        wnr, wnf = to_coords(wn_sq)
                        dr_n, df_n = abs(bkr - wnr), abs(bkf - wnf)
                        if (dr_n == 1 and df_n == 2) or (dr_n == 2 and df_n == 1):
                            check_by_knight = True
                        
                        # Pawn check?
                        if bkr == wpr + 1 and abs(bkf - wpf) == 1:
                            check_by_pawn = True
                        
                        if check_by_knight and check_by_pawn:
                            double_check_mates += 1
                        elif check_by_knight:
                            knight_mates += 1
                        elif check_by_pawn:
                            pawn_mates += 1

    total_mates = knight_mates + pawn_mates + double_check_mates
    
    # Print the final results in an equation format
    print(f"Checkmates by Knight alone: {knight_mates}")
    print(f"Checkmates by Pawn alone: {pawn_mates}")
    print(f"Checkmates by Double Check (Knight and Pawn): {double_check_mates}")
    print(f"Total Checkmate Positions = {knight_mates} + {pawn_mates} + {double_check_mates} = {total_mates}")

if __name__ == '__main__':
    solve_chess_puzzle()
