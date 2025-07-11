import itertools

def find_checkmate_positions():
    """
    Calculates the number of legal checkmate positions with White (King, Knight, Pawn)
    vs Black (King).

    This function iterates through all possible unique placements of the four pieces on an
    8x8 chessboard and checks each position for legality and checkmate conditions.
    """

    # Helper function to convert a square index (0-63) to board coordinates (row, col)
    # Board is indexed from top-left (0,0) = a8 to bottom-right (7,7) = h1.
    def to_coords(sq):
        return sq // 8, sq % 8

    # Helper function to get all squares attacked by a king at a given coordinate.
    def get_king_attacks(r, c):
        attacks = set()
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0:
                    continue
                nr, nc = r + dr, c + dc
                if 0 <= nr < 8 and 0 <= nc < 8:
                    attacks.add(nr * 8 + nc)
        return attacks

    # Helper function to get all squares attacked by a knight.
    def get_knight_attacks(r, c):
        attacks = set()
        moves = [(-2, -1), (-2, 1), (-1, -2), (-1, 2),
                 (1, -2), (1, 2), (2, -1), (2, 1)]
        for dr, dc in moves:
            nr, nc = r + dr, c + dc
            if 0 <= nr < 8 and 0 <= nc < 8:
                attacks.add(nr * 8 + nc)
        return attacks

    # Helper function to get squares attacked by a white pawn.
    # White pawns move "up" the board (row index decreases).
    def get_pawn_attacks(r, c):
        attacks = set()
        nr = r - 1
        if 0 <= nr < 8:
            if c > 0:  # Attack diagonally left
                attacks.add(nr * 8 + (c - 1))
            if c < 7:  # Attack diagonally right
                attacks.add(nr * 8 + (c + 1))
        return attacks

    checkmate_count = 0
    squares = range(64)
    
    # Using four nested loops is more efficient here than itertools.permutations
    # as we can filter illegal positions in the outer loops to prune the search space.
    
    # Loop over pawn positions first to apply the rank filter early.
    for wp_sq in squares:
        wp_r, wp_c = to_coords(wp_sq)
        # Legality: Pawns cannot be on the 1st (row 7) or 8th (row 0) rank.
        if wp_r == 0 or wp_r == 7:
            continue
        wp_attack_set = get_pawn_attacks(wp_r, wp_c)

        for wk_sq in squares:
            if wk_sq == wp_sq: continue
            wk_r, wk_c = to_coords(wk_sq)
            wk_attack_set = get_king_attacks(wk_r, wk_c)
            
            for bk_sq in squares:
                if bk_sq == wp_sq or bk_sq == wk_sq: continue
                bk_r, bk_c = to_coords(bk_sq)
                # Legality: Kings cannot be adjacent.
                if max(abs(wk_r - bk_r), abs(wk_c - bk_c)) <= 1:
                    continue
                
                # The black king's possible escape squares.
                bk_escape_routes = get_king_attacks(bk_r, bk_c)
                
                for wn_sq in squares:
                    if wn_sq in {wp_sq, wk_sq, bk_sq}: continue
                    
                    # ---- At this point, a unique, legal placement is defined ----
                    wn_r, wn_c = to_coords(wn_sq)
                    wn_attack_set = get_knight_attacks(wn_r, wn_c)
                    
                    all_white_attacks = wk_attack_set | wn_attack_set | wp_attack_set

                    # Check 1: Is the Black King in check?
                    if bk_sq not in all_white_attacks:
                        continue
                    
                    # Check 2: Does the Black King have any legal moves?
                    has_legal_move = False
                    white_piece_squares = {wk_sq, wn_sq, wp_sq}
                    for move_sq in bk_escape_routes:
                        # A move is legal if the destination is not occupied AND not attacked.
                        if move_sq not in white_piece_squares and move_sq not in all_white_attacks:
                            has_legal_move = True
                            break  # Found a safe square, so not a checkmate.
                            
                    # If no legal moves were found, it's a checkmate.
                    if not has_legal_move:
                        checkmate_count += 1
                        
    # The problem asks to output the number in the final equation.
    print(f"Number of checkmate positions = {checkmate_count}")
    return checkmate_count

# Run the function and capture the final answer.
# Note: This script is computationally intensive and may take a few minutes to run.
final_answer = find_checkmate_positions()
print(f"<<<{final_answer}>>>")