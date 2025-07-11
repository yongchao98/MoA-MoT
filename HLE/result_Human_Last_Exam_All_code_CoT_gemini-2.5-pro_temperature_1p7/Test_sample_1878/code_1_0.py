import itertools

def solve_chess_checkmates():
    """
    This script calculates the number of legal checkmate positions in a chess endgame
    where White has a King, a Knight, and a Pawn, and Black has only a King.

    The plan is to check every possible placement of these four pieces on the board
    and test each position for legality and checkmate conditions.
    """

    # We define a board coordinate system where a square is represented by (row, col).
    # (0, 0) corresponds to square a8, and (7, 7) corresponds to h1.
    # White pawns move from higher row index to lower (e.g., from row 6 to 5).

    def get_king_attacks(r, c):
        """Returns a set of squares {(r, c), ...} attacked by a king."""
        attacks = set()
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0:
                    continue
                nr, nc = r + dr, c + dc
                if 0 <= nr <= 7 and 0 <= nc <= 7:
                    attacks.add((nr, nc))
        return attacks

    def get_knight_attacks(r, c):
        """Returns a set of squares attacked by a knight."""
        attacks = set()
        moves = [(1, 2), (1, -2), (-1, 2), (-1, -2),
                 (2, 1), (2, -1), (-2, 1), (-2, -1)]
        for dr, dc in moves:
            nr, nc = r + dr, c + dc
            if 0 <= nr <= 7 and 0 <= nc <= 7:
                attacks.add((nr, nc))
        return attacks

    def get_pawn_attacks(r, c):
        """Returns a set of squares attacked by a white pawn."""
        attacks = set()
        nr = r - 1  # Pawns move "up" the board (row index decreases)
        if 0 <= nr <= 7:
            for dc in [-1, 1]:
                nc = c + dc
                if 0 <= nc <= 7:
                    attacks.add((nr, nc))
        return attacks

    mate_count = 0
    # Generate all unique placements of the 4 pieces on the 64 squares.
    # This generates over 230 million permutations to check.
    piece_placements = itertools.permutations(range(64), 4)

    for wk_sq, wn_sq, wp_sq, bk_sq in piece_placements:
        # Convert square indices (0-63) to (row, col) coordinates.
        wk_r, wk_c = divmod(wk_sq, 8)
        wn_r, wn_c = divmod(wn_sq, 8)
        wp_r, wp_c = divmod(wp_sq, 8)
        bk_r, bk_c = divmod(bk_sq, 8)

        # --- Legality Filter 1: Pawn on an invalid rank ---
        # A white pawn cannot be on rank 1 (row 7) or rank 8 (row 0).
        if not (1 <= wp_r <= 6):
            continue

        # --- Legality Filter 2: Kings are adjacent ---
        # The distance between kings must be > 1 in any direction (Chebyshev distance).
        if max(abs(wk_r - bk_r), abs(wk_c - bk_c)) <= 1:
            continue

        # --- Checkmate Analysis ---

        # 1. Determine attacks by the Knight and Pawn.
        knight_attacks = get_knight_attacks(wn_r, wn_c)
        pawn_attacks = get_pawn_attacks(wp_r, wp_c)
        
        # 2. Check if the Black king is currently in check.
        black_king_pos = (bk_r, bk_c)
        is_in_check = black_king_pos in knight_attacks or black_king_pos in pawn_attacks
        
        if not is_in_check:
            continue
            
        # 3. Check if the Black king has any legal moves.
        # First, find all squares controlled by *all* white pieces.
        king_attacks = get_king_attacks(wk_r, wk_c)
        all_white_attacks = king_attacks | knight_attacks | pawn_attacks
        
        # Next, find all potential escape squares for the Black king.
        black_king_potential_moves = get_king_attacks(bk_r, bk_c)
        
        # Legal moves are potential moves to squares not under attack.
        legal_black_king_moves = black_king_potential_moves - all_white_attacks
        
        # If the king is in check and has no legal moves, it's checkmate.
        if len(legal_black_king_moves) == 0:
            mate_count += 1
            
    print(f"The number of possible checkmate positions is: {mate_count}")

# Note: This script performs an exhaustive search and may take several minutes to complete.
if __name__ == '__main__':
    solve_chess_checkmates()