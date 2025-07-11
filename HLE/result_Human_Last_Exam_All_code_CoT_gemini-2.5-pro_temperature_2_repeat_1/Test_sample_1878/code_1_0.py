import itertools
import time

def to_coords(sq):
    """Converts a square index (0-63) to (row, col) coordinates."""
    return sq // 8, sq % 8

def to_sq(r, c):
    """Converts (row, col) coordinates to a square index (0-63)."""
    return r * 8 + c

def get_king_attacks(sq, king_attack_map):
    """Returns a set of squares a king attacks from a given square."""
    return king_attack_map.get(sq, set())

def get_knight_attacks(sq, knight_attack_map):
    """Returns a set of squares a knight attacks from a given square."""
    return knight_attack_map.get(sq, set())

def get_pawn_attacks(sq, pawn_attack_map):
    """Returns a set of squares a white pawn attacks from a given square."""
    return pawn_attack_map.get(sq, set())

def precompute_attacks():
    """Precomputes attack sets for all pieces on all squares to speed up checks."""
    king_attack_map = {}
    knight_attack_map = {}
    pawn_attack_map = {}
    
    for r in range(8):
        for c in range(8):
            sq = to_sq(r, c)
            
            # King attacks
            k_attacks = set()
            for dr in [-1, 0, 1]:
                for dc in [-1, 0, 1]:
                    if dr == 0 and dc == 0: continue
                    nr, nc = r + dr, c + dc
                    if 0 <= nr < 8 and 0 <= nc < 8:
                        k_attacks.add(to_sq(nr, nc))
            king_attack_map[sq] = k_attacks

            # Knight attacks
            n_attacks = set()
            for dr, dc in [(1, 2), (1, -2), (-1, 2), (-1, -2), (2, 1), (2, -1), (-2, 1), (-2, -1)]:
                nr, nc = r + dr, c + dc
                if 0 <= nr < 8 and 0 <= nc < 8:
                    n_attacks.add(to_sq(nr, nc))
            knight_attack_map[sq] = n_attacks

            # Pawn attacks (White moves from row 0 to 7)
            p_attacks = set()
            if r < 7: # Pawns cannot attack from the 8th rank
                for dc in [-1, 1]:
                    nr, nc = r + 1, c + dc
                    if 0 <= nr < 8 and 0 <= nc < 8:
                        p_attacks.add(to_sq(nr, nc))
            pawn_attack_map[sq] = p_attacks
            
    return king_attack_map, knight_attack_map, pawn_attack_map

def solve_chess_puzzle():
    """
    Finds the number of legal checkmate positions with WK, WN, WP vs BK.
    """
    print("Pre-computing attack patterns...")
    king_attack_map, knight_attack_map, pawn_attack_map = precompute_attacks()
    
    count = 0
    squares = range(64)
    total_permutations = 64 * 63 * 62 * 61
    
    print(f"Starting search over {total_permutations:,} positions. This may take a few minutes...")
    start_time = time.time()
    
    # Iterate through all permutations of 4 unique squares for the 4 pieces
    for i, piece_combo in enumerate(itertools.permutations(squares, 4)):
        wk_sq, wn_sq, wp_sq, bk_sq = piece_combo
        
        # Performance update
        if i > 0 and i % 2000000 == 0:
            elapsed = time.time() - start_time
            print(f"  ...processed {i / total_permutations:.0%} ({i:,} positions) in {elapsed:.1f}s")

        # --- Legality Checks ---
        wp_r, _ = to_coords(wp_sq)
        # 1. White pawn cannot be on rank 1 (row 0) or rank 8 (row 7)
        if wp_r == 0 or wp_r == 7:
            continue
            
        # 2. Kings cannot be on adjacent squares
        if bk_sq in get_king_attacks(wk_sq, king_attack_map):
            continue

        # --- Checkmate Checks ---
        # 1. Is the Black King in check?
        is_checked_by_knight = bk_sq in get_knight_attacks(wn_sq, knight_attack_map)
        is_checked_by_pawn = bk_sq in get_pawn_attacks(wp_sq, pawn_attack_map)
        
        if not (is_checked_by_knight or is_checked_by_pawn):
            continue

        # 2. Can the Black King escape?
        # A move is legal if the destination square is not occupied and not attacked.
        white_pieces = {wk_sq, wn_sq, wp_sq}
        white_attacks = (get_king_attacks(wk_sq, king_attack_map) |
                         get_knight_attacks(wn_sq, knight_attack_map) |
                         get_pawn_attacks(wp_sq, pawn_attack_map))

        has_escape_move = False
        for escape_sq in get_king_attacks(bk_sq, king_attack_map):
            if escape_sq not in white_pieces and escape_sq not in white_attacks:
                has_escape_move = True
                break
        
        if has_escape_move:
            continue
            
        # If we reach here, Black is in check and has no legal moves. It's checkmate.
        count += 1
        
    end_time = time.time()
    print(f"\nSearch complete in {end_time - start_time:.2f} seconds.")

    # The problem asks to output the numbers in the "final equation".
    # This is interpreted as printing the components of the result.
    print("\nResult:")
    print(f"Total possible checkmate positions found: {count}")

if __name__ == '__main__':
    solve_chess_puzzle()

<<<4156>>>