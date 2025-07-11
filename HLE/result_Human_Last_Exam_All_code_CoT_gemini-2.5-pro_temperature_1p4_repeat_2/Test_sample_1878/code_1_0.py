import itertools

def solve_chess_checkmates():
    """
    Calculates the number of legal checkmate positions where White has a king,
    pawn, and knight, and Black has only a king.
    """
    count = 0
    squares = range(64)

    # Pre-calculate attack patterns for each piece from each square
    # Board representation: a1=0, b1=1, ..., h1=7, a2=8, ..., h8=63
    # A square's coordinate (row, col) is (sq // 8, sq % 8)
    # White pawns move from lower row index to higher row index.
    king_moves = {}
    knight_moves = {}
    pawn_attacks = {}

    for sq in squares:
        r, c = sq // 8, sq % 8
        
        # King moves
        k_moves = set()
        for dr, dc in itertools.product([-1, 0, 1], repeat=2):
            if dr == 0 and dc == 0: continue
            nr, nc = r + dr, c + dc
            if 0 <= nr < 8 and 0 <= nc < 8:
                k_moves.add(nr * 8 + nc)
        king_moves[sq] = k_moves
        
        # Knight moves
        n_moves = set()
        for dr, dc in [(1, 2), (1, -2), (-1, 2), (-1, -2), (2, 1), (2, -1), (-2, 1), (-2, -1)]:
            nr, nc = r + dr, c + dc
            if 0 <= nr < 8 and 0 <= nc < 8:
                n_moves.add(nr * 8 + nc)
        knight_moves[sq] = n_moves
        
        # White Pawn attacks
        p_attacks = set()
        # A pawn at (r, c) attacks (r+1, c-1) and (r+1, c+1)
        nr = r + 1
        if 0 <= nr < 8:
            for dc in [-1, 1]:
                nc = c + dc
                if 0 <= nc < 8:
                    p_attacks.add(nr * 8 + nc)
        pawn_attacks[sq] = p_attacks

    # Iterate through all possible unique placements of the 4 pieces
    for wk_sq, bk_sq, wn_sq, wp_sq in itertools.permutations(squares, 4):
        
        # --- Legality Checks ---
        
        # 1. White pawn cannot be on the 1st or 8th rank (row 0 or 7)
        wp_r = wp_sq // 8
        if wp_r == 0 or wp_r == 7:
            continue
            
        # 2. Kings cannot be on adjacent squares
        wk_r, wk_c = wk_sq // 8, wk_sq % 8
        bk_r, bk_c = bk_sq // 8, bk_sq % 8
        if max(abs(wk_r - bk_r), abs(wk_c - bk_c)) <= 1:
            continue

        # --- Checkmate Analysis ---
        white_piece_squares = {wk_sq, wn_sq, wp_sq}
        
        # Determine all squares attacked by White
        attacked_by_white = king_moves[wk_sq] | knight_moves[wn_sq] | pawn_attacks[wp_sq]

        # Condition 1: The Black King must be in check.
        if bk_sq not in attacked_by_white:
            continue

        # Condition 2: The Black King must have no legal moves.
        has_safe_escape = False
        for escape_sq in king_moves[bk_sq]:
            # A move is legal if the destination is not occupied by a friendly (White) piece
            # AND not under attack.
            if escape_sq not in white_piece_squares and escape_sq not in attacked_by_white:
                has_safe_escape = True
                break
        
        # If the king is in check and has no safe escape, it's checkmate.
        if not has_safe_escape:
            count += 1
            
    print(count)

if __name__ == '__main__':
    solve_chess_checkmates()