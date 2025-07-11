import itertools
import sys

def solve_chess_problem():
    """
    This script calculates the number of legal checkmate positions with 
    White (King, Pawn, Knight) vs Black (King).
    
    A position is 'legal' if it can be reached from a position where the 
    side not to move (Black) was not in check.
    """
    
    # SETUP: Pre-calculate move patterns for all 64 squares
    # We represent the board with (row, col) tuples, where (0,0) is a1 and (7,7) is h8.
    
    def _get_raw_knight_moves(r, c):
        moves = []
        deltas = [(-2, -1), (-2, 1), (-1, -2), (-1, 2), (1, -2), (1, 2), (2, -1), (2, 1)]
        for dr, dc in deltas:
            nr, nc = r + dr, c + dc
            if 0 <= nr <= 7 and 0 <= nc <= 7:
                moves.append((nr, nc))
        return moves

    def _get_raw_king_moves(r, c):
        moves = []
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0: continue
                nr, nc = r + dr, c + dc
                if 0 <= nr <= 7 and 0 <= nc <= 7:
                    moves.append((nr, nc))
        return moves

    ALL_SQUARES = [(r, c) for r in range(8) for c in range(8)]
    KNIGHT_MOVES = {pos: _get_raw_knight_moves(*pos) for pos in ALL_SQUARES}
    KING_MOVES = {pos: _get_raw_king_moves(*pos) for pos in ALL_SQUARES}

    def get_pawn_attacks(pawn_pos):
        r, c = pawn_pos
        attacks = []
        if r < 7:  # A pawn on the 8th rank would be promoted.
            if c > 0: attacks.append((r + 1, c - 1))
            if c < 7: attacks.append((r + 1, c + 1))
        return attacks

    # --- HELPER FUNCTIONS ---
    
    def is_square_attacked_by_white(target_sq, wk, wp, wn):
        """Checks if a square is geometrically attacked by any white piece."""
        if wn and target_sq in KNIGHT_MOVES.get(wn, []): return True
        if wp and target_sq in get_pawn_attacks(wp): return True
        if wk and target_sq in KING_MOVES.get(wk, []): return True
        return False

    def find_white_attackers(target_sq, wk, wp, wn):
        """Finds which white pieces are attacking a square."""
        attackers = []
        if wn and target_sq in KNIGHT_MOVES.get(wn, []): attackers.append(wn)
        if wp and target_sq in get_pawn_attacks(wp): attackers.append(wp)
        if wk and target_sq in KING_MOVES.get(wk, []): attackers.append(wk)
        return attackers

    # --- CORE LOGIC FUNCTIONS ---

    def is_checkmate(wk, bk, wp, wn):
        """Determines if the black king is checkmated."""
        attackers = find_white_attackers(bk, wk, wp, wn)
        if not attackers:
            return False

        # Check for any legal move for the black king.
        for move in KING_MOVES[bk]:
            # A king cannot move adjacent to the other king.
            if max(abs(move[0] - wk[0]), abs(move[1] - wk[1])) <= 1:
                continue

            # Can the king move to an empty, unattacked square?
            if move not in [wk, wp, wn]:
                if not is_square_attacked_by_white(move, wk, wp, wn):
                    return False  # Found a safe empty square.

            # Can the king capture an attacking piece?
            elif move in attackers and len(attackers) == 1:
                # To capture, the destination must be safe from other pieces.
                temp_wk, temp_wp, temp_wn = wk, wp, wn
                if move == wk: temp_wk = None
                elif move == wp: temp_wp = None
                elif move == wn: temp_wn = None
                
                if not is_square_attacked_by_white(move, temp_wk, temp_wp, temp_wn):
                    return False  # King can safely capture the checker.
        
        return True # No legal moves found, so it is checkmate.

    def is_legal_by_retroanalysis(wk, bk, wp, wn):
        """Checks if the position could have arisen from a legal previous position."""
        white_pieces = {wk, wp, wn}

        # Try to un-move the Knight
        for prev_pos in KNIGHT_MOVES[wn]:
            if prev_pos not in white_pieces and prev_pos != bk:
                if not is_square_attacked_by_white(bk, wk, wp, prev_pos): return True

        # Try to un-move the King
        for prev_pos in KING_MOVES[wk]:
            if prev_pos not in white_pieces and prev_pos != bk:
                if max(abs(prev_pos[0] - bk[0]), abs(prev_pos[1] - bk[1])) > 1:
                    if not is_square_attacked_by_white(bk, prev_pos, wp, wn): return True

        # Try to un-move the Pawn
        wp_r, wp_c = wp
        # One-step move
        if wp_r > 0:
            prev_pos = (wp_r - 1, wp_c)
            if prev_pos not in white_pieces and prev_pos != bk:
                if not is_square_attacked_by_white(bk, wk, prev_pos, wn): return True
        # Two-step move (from rank 2 to 4)
        if wp_r == 3:
            prev_pos = (1, wp_c)
            path_sq = (2, wp_c)
            if prev_pos not in white_pieces and prev_pos != bk and path_sq not in white_pieces and path_sq != bk:
                if not is_square_attacked_by_white(bk, wk, prev_pos, wn): return True

        return False

    # --- MAIN EXECUTION ---
    
    total_permutations = 64 * 63 * 62 * 61
    valid_setup_count = 0
    checkmate_count = 0
    legal_checkmate_count = 0
    
    all_squares_1d = range(64)
    
    print("Starting analysis... This may take several minutes.")
    
    for i, p in enumerate(itertools.permutations(all_squares_1d, 4)):
        if i % 1000000 == 0 and i > 0:
            print(f"  ...processed {i:,} of {total_permutations:,} permutations...")

        wk_pos = (p[0] // 8, p[0] % 8)
        bk_pos = (p[1] // 8, p[1] % 8)
        wp_pos = (p[2] // 8, p[2] % 8)
        wn_pos = (p[3] // 8, p[3] % 8)

        # Constraint 1: Basic setup validity (Pawn on ranks 2-7, Kings not adjacent)
        if not (1 <= wp_pos[0] <= 6) or max(abs(wk_pos[0] - bk_pos[0]), abs(wk_pos[1] - bk_pos[1])) <= 1:
            continue
        valid_setup_count += 1

        # Constraint 2: Must be a checkmate position
        if not is_checkmate(wk_pos, bk_pos, wp_pos, wn_pos):
            continue
        checkmate_count += 1
        
        # Constraint 3: Must be reachable from a legal position (retro-analysis)
        if not is_legal_by_retroanalysis(wk_pos, bk_pos, wp_pos, wn_pos):
            continue
        legal_checkmate_count += 1
    
    print("\nAnalysis complete. Here are the results:\n")
    print("-" * 60)
    print("Step 1: Total possible placements of the 4 pieces.")
    print(f"Total permutations = 64 * 63 * 62 * 61 = {total_permutations:,}")
    print("\nStep 2: Filter for basic board rules (pawn on ranks 2-7, kings not adjacent).")
    print(f"Number of valid setups found = {valid_setup_count:,}")
    print("\nStep 3: From those, filter for positions that are checkmate.")
    print(f"Number of checkmate positions found = {checkmate_count:,}")
    print("\nStep 4: From those, filter for positions that are 'provably legal' via retro-analysis.")
    print(f"Number of legal checkmate positions found = {legal_checkmate_count:,}")
    print("-" * 60)
    print("The final number of constructible legal checkmates is the result of the final step.")
    print(f"Final Equation: The number of legal positions is {legal_checkmate_count}")
    
    return legal_checkmate_count

# Run the calculation and print the final answer in the requested format.
final_answer = solve_chess_problem()
print(f"\n<<<{final_answer}>>>")
