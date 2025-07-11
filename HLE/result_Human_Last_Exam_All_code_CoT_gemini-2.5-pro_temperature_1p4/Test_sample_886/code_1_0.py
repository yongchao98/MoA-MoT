import collections

def solve():
    """
    Calculates the number of distinct checkmate positions for hybrid chess pieces.

    A hybrid piece is a combination of two standard pieces. This script checks
    the most powerful combinations:
    - Chancellor (Rook + Knight)
    - Princess (Bishop + Knight)
    - Amazon (Queen + Knight, which is Rook + Bishop + Knight)

    A checkmate position is one where the king is under attack and all its
    escape squares are also attacked by the hybrid piece. An adjacent piece
    is not a mate, as the king could capture it.

    Distinct positions are counted by normalizing each found mate through
    board symmetries (rotations and reflections) and counting the unique
    canonical results.
    """
    
    # --- Precompute piece attacks and king moves for all 64 squares ---
    ROOK_ATTACKS = collections.defaultdict(set)
    BISHOP_ATTACKS = collections.defaultdict(set)
    KNIGHT_ATTACKS = collections.defaultdict(set)
    KING_MOVES = collections.defaultdict(set)

    for r in range(8):
        for c in range(8):
            # King
            for dr in [-1, 0, 1]:
                for dc in [-1, 0, 1]:
                    if dr == 0 and dc == 0: continue
                    nr, nc = r + dr, c + dc
                    if 0 <= nr < 8 and 0 <= nc < 8:
                        KING_MOVES[(r, c)].add((nr, nc))
            # Rook
            for i in range(8):
                if i != c: ROOK_ATTACKS[(r, c)].add((r, i))
                if i != r: ROOK_ATTACKS[(r, c)].add((i, c))
            # Bishop
            for i in range(1, 8):
                if 0 <= r + i < 8 and 0 <= c + i < 8: BISHOP_ATTACKS[(r, c)].add((r + i, c + i))
                if 0 <= r - i < 8 and 0 <= c + i < 8: BISHOP_ATTACKS[(r, c)].add((r - i, c + i))
                if 0 <= r + i < 8 and 0 <= c - i < 8: BISHOP_ATTACKS[(r, c)].add((r + i, c - i))
                if 0 <= r - i < 8 and 0 <= c - i < 8: BISHOP_ATTACKS[(r, c)].add((r - i, c - i))
            # Knight
            for dr, dc in [(1, 2), (1, -2), (-1, 2), (-1, -2), (2, 1), (2, -1), (-2, 1), (-2, -1)]:
                nr, nc = r + dr, c + dc
                if 0 <= nr < 8 and 0 <= nc < 8:
                    KNIGHT_ATTACKS[(r, c)].add((nr, nc))

    def get_hybrid_attacks(piece_type, pos):
        attacks = set()
        if 'R' in piece_type: attacks.update(ROOK_ATTACKS[pos])
        if 'B' in piece_type: attacks.update(BISHOP_ATTACKS[pos])
        if 'N' in piece_type: attacks.update(KNIGHT_ATTACKS[pos])
        return attacks

    def get_canonical(pos):
        k_pos, p_pos = tuple(sorted(pos))
        symmetries = set()
        ck, cp = k_pos, p_pos
        for _ in range(4): # 4 rotations
            symmetries.add(tuple(sorted(((ck[0], ck[1]), (cp[0], cp[1]))))) # Add original
            symmetries.add(tuple(sorted(((ck[0], 7 - ck[1]), (cp[0], 7 - cp[1]))))) # Add reflection
            ck = (ck[1], 7 - ck[0]) # Rotate 90 degrees
            cp = (cp[1], 7 - cp[0])
        return min(symmetries)

    hybrid_pieces = {
        "Chancellor (Rook + Knight)": "RN",
        "Princess (Bishop + Knight)": "BN",
        "Amazon (Queen + Knight)": "RBN" # Queen = Rook + Bishop
    }

    total_distinct_mates = 0
    final_equations = []

    for name, piece_code in hybrid_pieces.items():
        canonical_mates = set()
        for kr in range(8):
            for kc in range(8):
                king_pos = (kr, kc)
                for pr in range(8):
                    for pc in range(8):
                        piece_pos = (pr, pc)
                        
                        if king_pos == piece_pos: continue
                        if piece_pos in KING_MOVES[king_pos]: continue

                        attacked_by_piece = get_hybrid_attacks(piece_code, piece_pos)
                        
                        if king_pos not in attacked_by_piece: continue
                            
                        can_escape = any(esc_pos not in attacked_by_piece for esc_pos in KING_MOVES[king_pos])
                        
                        if not can_escape:
                            canonical_mates.add(get_canonical((king_pos, piece_pos)))
        
        num_distinct = len(canonical_mates)
        final_equations.append(f"{name}: {num_distinct}")
        total_distinct_mates += num_distinct

    print("Number of distinct checkmate positions for each hybrid piece:")
    for eq in final_equations:
        print(eq)
    
    print("\nTotal distinct checkmate positions:")
    nums = [int(s.split(': ')[1]) for s in final_equations]
    sum_str = " + ".join(map(str, nums))
    print(f"{sum_str} = {total_distinct_mates}")


if __name__ == '__main__':
    solve()
    print("\n<<<4>>>")