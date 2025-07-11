import itertools

def solve_chess_puzzle():
    """
    Calculates the number of legal checkmate positions with K, N, P vs K.
    """
    
    # --- Pre-computation of piece moves for efficiency ---
    KING_MOVES = [[] for _ in range(64)]
    KNIGHT_MOVES = [[] for _ in range(64)]
    # PAWN_ATTACKS[sq] lists the squares a white pawn on sq attacks.
    PAWN_ATTACKS = [[] for _ in range(64)]

    def coord(sq):
        return sq // 8, sq % 8

    def sq_from_coord(r, c):
        return r * 8 + c

    for r_start in range(8):
        for c_start in range(8):
            sq = sq_from_coord(r_start, c_start)
            
            # Pre-compute King moves
            for r_delta in [-1, 0, 1]:
                for c_delta in [-1, 0, 1]:
                    if r_delta == 0 and c_delta == 0: continue
                    r_end, c_end = r_start + r_delta, c_start + c_delta
                    if 0 <= r_end < 8 and 0 <= c_end < 8:
                        KING_MOVES[sq].append(sq_from_coord(r_end, c_end))

            # Pre-compute Knight moves
            for r_delta, c_delta in [(1,2), (1,-2), (-1,2), (-1,-2), (2,1), (2,-1), (-2,1), (-2,-1)]:
                r_end, c_end = r_start + r_delta, c_start + c_delta
                if 0 <= r_end < 8 and 0 <= c_end < 8:
                    KNIGHT_MOVES[sq].append(sq_from_coord(r_end, c_end))
            
            # Pre-compute White Pawn attacks
            if r_start < 7: # A pawn on the 8th rank cannot attack
                if c_start > 0: PAWN_ATTACKS[sq].append(sq_from_coord(r_start + 1, c_start - 1))
                if c_start < 7: PAWN_ATTACKS[sq].append(sq_from_coord(r_start + 1, c_start + 1))

    # --- Helper Functions ---
    def is_attacked_by_white(target_sq, wk_sq, wn_sq, wp_sq):
        """Checks if a square is attacked by any of the specified white pieces."""
        if wk_sq is not None and target_sq in KING_MOVES[wk_sq]: return True
        if wn_sq is not None and target_sq in KNIGHT_MOVES[wn_sq]: return True
        if wp_sq is not None and target_sq in PAWN_ATTACKS[wp_sq]: return True
        return False

    def is_legal_mate_position(wk_sq, bk_sq, wn_sq, wp_sq):
        """
        Checks if the position could have arisen from a legal previous move.
        This is done by 'un-moving' each white piece and checking if the black
        king was NOT in check in the resulting position.
        """
        # Un-move the White Knight
        occupied_without_n = {wk_sq, bk_sq, wp_sq}
        for from_sq in KNIGHT_MOVES[wn_sq]:
            if from_sq not in occupied_without_n:
                if not is_attacked_by_white(bk_sq, wk_sq, from_sq, wp_sq): return True

        # Un-move the White Pawn
        occupied_without_p = {wk_sq, bk_sq, wn_sq}
        wp_row, wp_col = coord(wp_sq)
        # Un-move 1 step
        if wp_row > 1:
            from_sq = sq_from_coord(wp_row - 1, wp_col)
            if from_sq not in occupied_without_p:
                if not is_attacked_by_white(bk_sq, wk_sq, wn_sq, from_sq): return True
        # Un-move 2 steps (pawn must have started on rank 2 and landed on rank 4)
        if wp_row == 3:
            from_sq = sq_from_coord(1, wp_col)
            mid_sq = sq_from_coord(2, wp_col)
            if from_sq not in occupied_without_p and mid_sq not in occupied_without_p:
                if not is_attacked_by_white(bk_sq, wk_sq, wn_sq, from_sq): return True

        # Un-move the White King
        occupied_without_k = {bk_sq, wn_sq, wp_sq}
        for from_sq in KING_MOVES[wk_sq]:
            # King cannot move out of check from the other king
            if from_sq not in KING_MOVES[bk_sq] and from_sq not in occupied_without_k:
                if not is_attacked_by_white(bk_sq, from_sq, wn_sq, wp_sq): return True
        
        return False

    # --- Main Logic ---
    count = 0
    all_squares = range(64)
    
    # Iterate through every possible placement of the 4 pieces on 64 squares
    for wk_sq, bk_sq, wn_sq, wp_sq in itertools.permutations(all_squares, 4):
        
        # 1. Filter based on basic placement rules
        if coord(wp_sq)[0] in (0, 7): continue  # Pawn can't be on 1st or 8th rank
        if wk_sq in KING_MOVES[bk_sq]: continue  # Kings can't be adjacent

        # 2. Filter based on checkmate conditions
        # 2a. Black king must be in check by exactly one piece
        attackers = []
        if bk_sq in KNIGHT_MOVES[wn_sq]: attackers.append(wn_sq)
        if bk_sq in PAWN_ATTACKS[wp_sq]: attackers.append(wp_sq)
        if len(attackers) != 1: continue
        checker_sq = attackers[0]

        # 2b. Black king must not have any escape squares
        has_escape = False
        occupied_squares = {wk_sq, wn_sq, wp_sq}
        for escape_sq in KING_MOVES[bk_sq]:
            if escape_sq in occupied_squares: continue
            if not is_attacked_by_white(escape_sq, wk_sq, wn_sq, wp_sq):
                has_escape = True
                break
        if has_escape: continue

        # 2c. The checking piece cannot be capturable by the black king
        can_capture = False
        if checker_sq in KING_MOVES[bk_sq]: # Is the checker adjacent?
            # Is the checker protected? Check if its square is attacked by other white pieces
            if checker_sq == wn_sq:
                if not is_attacked_by_white(checker_sq, wk_sq, None, wp_sq):
                    can_capture = True
            elif checker_sq == wp_sq:
                if not is_attacked_by_white(checker_sq, wk_sq, wn_sq, None):
                    can_capture = True
        if can_capture: continue
        
        # 3. Filter based on position legality (retrograde analysis)
        if is_legal_mate_position(wk_sq, bk_sq, wn_sq, wp_sq):
            count += 1
            
    print(f"{count}")

solve_chess_puzzle()
<<<468>>>