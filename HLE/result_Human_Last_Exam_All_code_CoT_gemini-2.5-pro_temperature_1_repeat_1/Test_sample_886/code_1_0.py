import itertools
from functools import lru_cache

def solve_chess_checkmates():
    """
    Calculates the number of distinct checkmate positions for various combined-move chess pieces.

    This script iterates through all possible board setups for a given white attacker,
    white king, and black king. For each setup, it checks if it constitutes a checkmate.
    To count only distinct positions, it normalizes each found mate position by checking
    all 8 board symmetries (rotations and reflections) and using the lexicographically
    smallest one as the canonical representation.
    """

    # Memoization to speed up repeated calculations for attack sets
    @lru_cache(maxsize=None)
    def is_valid(r, c):
        return 0 <= r < 8 and 0 <= c < 8

    @lru_cache(maxsize=64)
    def get_rook_attacks(pos):
        r, c = pos
        attacks = set()
        for dr, dc in [(1, 0), (-1, 0), (0, 1), (0, -1)]:
            for i in range(1, 8):
                nr, nc = r + i * dr, c + i * dc
                if is_valid(nr, nc):
                    attacks.add((nr, nc))
                else:
                    break
        return attacks

    @lru_cache(maxsize=64)
    def get_bishop_attacks(pos):
        r, c = pos
        attacks = set()
        for dr, dc in [(1, 1), (1, -1), (-1, 1), (-1, -1)]:
            for i in range(1, 8):
                nr, nc = r + i * dr, c + i * dc
                if is_valid(nr, nc):
                    attacks.add((nr, nc))
                else:
                    break
        return attacks

    @lru_cache(maxsize=64)
    def get_knight_attacks(pos):
        r, c = pos
        attacks = set()
        for dr, dc in [(1, 2), (1, -2), (-1, 2), (-1, -2), (2, 1), (2, -1), (-2, 1), (-2, -1)]:
            nr, nc = r + dr, c + dc
            if is_valid(nr, nc):
                attacks.add((nr, nc))
        return attacks

    @lru_cache(maxsize=64)
    def get_king_attacks(pos):
        r, c = pos
        attacks = set()
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0:
                    continue
                nr, nc = r + dr, c + dc
                if is_valid(nr, nc):
                    attacks.add((nr, nc))
        return attacks

    # Define the unique, mate-capable pieces based on combinations.
    # The moves are defined by a tuple of base attack functions.
    PIECE_DEFINITIONS = {
        "Rook": (get_rook_attacks,),
        "Queen (B+R)": (get_rook_attacks, get_bishop_attacks),
        "Archbishop (B+N)": (get_bishop_attacks, get_knight_attacks),
        "Chancellor (R+N)": (get_rook_attacks, get_knight_attacks),
        "KnightKing (K+N)": (get_king_attacks, get_knight_attacks),
        "BishopKing (K+B)": (get_king_attacks, get_bishop_attacks),
        "RookKing (K+R)": (get_king_attacks, get_rook_attacks),
        "Amazon (Q+N)": (get_rook_attacks, get_bishop_attacks, get_knight_attacks),
        "QueenKing (K+Q)": (get_king_attacks, get_rook_attacks, get_bishop_attacks),
    }

    # Helper for canonical representation
    def get_canonical(p1, p2, p3):
        positions = [p1, p2, p3]
        symmetries = set()
        for _ in range(4):  # 4 rotations
            # Original and horizontal flip
            symmetries.add(tuple(sorted(positions)))
            flipped_positions = [(r, 7 - c) for r, c in positions]
            symmetries.add(tuple(sorted(flipped_positions)))
            # Rotate for next iteration
            positions = [(c, 7 - r) for r, c in positions]
        return min(symmetries)

    all_coords = list(itertools.product(range(8), repeat=2))
    piece_results = {}
    total_mates = 0

    print("Calculating distinct checkmate positions for each piece type...")
    print("This may take several minutes.")

    for name, attack_funcs in PIECE_DEFINITIONS.items():
        canonical_mates = set()
        # Iterate through all permutations of 3 distinct squares for the 3 pieces
        for wa_pos, wk_pos, bk_pos in itertools.permutations(all_coords, 3):

            # Rule: Kings cannot be on adjacent squares
            if abs(wk_pos[0] - bk_pos[0]) <= 1 and abs(wk_pos[1] - bk_pos[1]) <= 1:
                continue

            # 1. Calculate all squares attacked by the white pieces
            attacker_attacks = set()
            for func in attack_funcs:
                attacker_attacks.update(func(wa_pos))
            
            wk_attacks = get_king_attacks(wk_pos)
            all_white_attacks = attacker_attacks | wk_attacks

            # 2. Check if the Black King is in check
            if bk_pos not in attacker_attacks:
                continue

            # 3. Check if the Black King has any legal moves to escape
            can_escape = False
            for bk_move in get_king_attacks(bk_pos):
                if bk_move not in all_white_attacks:
                    can_escape = True
                    break
            
            # 4. If in check and cannot escape, it's checkmate
            if not can_escape:
                canon_form = get_canonical(wa_pos, wk_pos, bk_pos)
                canonical_mates.add(canon_form)

        count = len(canonical_mates)
        piece_results[name] = count
        total_mates += count
        print(f"  - Found {count} distinct positions for: {name}")

    print("\n--- Final Results ---")
    
    # Sort results for consistent output
    sorted_pieces = sorted(piece_results.keys())
    
    sum_parts = []
    for name in sorted_pieces:
        count = piece_results[name]
        print(f"{name}: {count}")
        sum_parts.append(str(count))

    equation = " + ".join(sum_parts)
    print(f"\nTotal Distinct Checkmate Positions = {equation} = {total_mates}")

solve_chess_checkmates()
<<<7441>>>