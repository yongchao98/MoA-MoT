import sys

def calculate_checkmates():
    """
    This script calculates the total number of legal checkmate positions where
    White has a King, Knight, and Pawn, and Black has only a King.
    """

    # Memoization caches to speed up move generation
    king_move_cache = {}
    knight_move_cache = {}
    pawn_attack_cache = {}

    def to_coords(sq):
        return sq // 8, sq % 8

    def get_king_moves(sq):
        if sq in king_move_cache:
            return king_move_cache[sq]
        moves = set()
        r, c = to_coords(sq)
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0:
                    continue
                nr, nc = r + dr, c + dc
                if 0 <= nr < 8 and 0 <= nc < 8:
                    moves.add(nr * 8 + nc)
        king_move_cache[sq] = moves
        return moves

    def get_knight_moves(sq):
        if sq in knight_move_cache:
            return knight_move_cache[sq]
        moves = set()
        r, c = to_coords(sq)
        deltas = [(-2, -1), (-2, 1), (-1, -2), (-1, 2),
                  (1, -2), (1, 2), (2, -1), (2, 1)]
        for dr, dc in deltas:
            nr, nc = r + dr, c + dc
            if 0 <= nr < 8 and 0 <= nc < 8:
                moves.add(nr * 8 + nc)
        knight_move_cache[sq] = moves
        return moves

    def get_pawn_attacks(sq):
        # Assuming white pawns move from low rank index to high rank index
        if sq in pawn_attack_cache:
            return pawn_attack_cache[sq]
        attacks = set()
        r, c = to_coords(sq)
        if r < 7:  # A pawn cannot attack from the last rank
            if c > 0: # Not on the 'a' file
                attacks.add((r + 1) * 8 + (c - 1))
            if c < 7: # Not on the 'h' file
                attacks.add((r + 1) * 8 + (c + 1))
        pawn_attack_cache[sq] = attacks
        return attacks

    def king_distance(sq1, sq2):
        r1, c1 = to_coords(sq1)
        r2, c2 = to_coords(sq2)
        return max(abs(r1 - r2), abs(c1 - c2))

    mate_count = 0
    # Main loop to iterate through all piece placements
    for bk_sq in range(64):
        # This loop is computationally intensive. Uncomment the following lines
        # to see progress printed to the console during execution.
        # if bk_sq % 4 == 0:
        #     print(f"Analyzing positions with Black King on square {bk_sq}/63...", file=sys.stderr)

        for wn_sq in range(64):
            if wn_sq == bk_sq: continue

            # Pawns can only be on ranks 2-7 (squares 8 through 55)
            for wp_sq in range(8, 56):
                if wp_sq == bk_sq or wp_sq == wn_sq: continue

                for wk_sq in range(64):
                    if wk_sq in {bk_sq, wn_sq, wp_sq}: continue

                    # --- Check for Mate ---
                    # 1. Legality: Kings cannot be adjacent.
                    if king_distance(wk_sq, bk_sq) <= 1:
                        continue

                    # 2. Check: Black King must be attacked by the Knight or Pawn.
                    wn_moves = get_knight_moves(wn_sq)
                    wp_attacks = get_pawn_attacks(wp_sq)
                    if bk_sq not in wn_moves and bk_sq not in wp_attacks:
                        continue

                    # 3. No Escape: Black King must have no legal moves.
                    has_legal_move = False
                    # Squares attacked by the white king are illegal for the black king to move to.
                    wk_attacks = get_king_moves(wk_sq)

                    for dest_sq in get_king_moves(bk_sq):
                        # A move is legal if the destination square is not attacked.
                        if dest_sq in wk_attacks: continue
                        if dest_sq in wn_moves: continue
                        if dest_sq in wp_attacks: continue

                        # If we reach here, we have found a safe escape square.
                        has_legal_move = True
                        break  # No need to check other moves

                    if not has_legal_move:
                        # This position is a valid checkmate.
                        mate_count += 1
                        
    # Final Output
    print(f"Total number of legal checkmate positions: {mate_count}")

calculate_checkmates()