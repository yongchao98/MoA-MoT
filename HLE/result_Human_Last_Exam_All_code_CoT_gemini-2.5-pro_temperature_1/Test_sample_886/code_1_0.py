import collections

def solve_chess_puzzle():
    """
    Calculates the total number of distinct checkmate positions achievable by
    a single piece that combines the moves of two standard chess pieces.
    """

    # Helper function to check if a square is attacked by a piece with combined moves.
    def is_attacked(attacker_pos, target_pos, piece_moves):
        if attacker_pos == target_pos:
            return False
        
        ra, ca = attacker_pos
        rt, ct = target_pos
        dr, dc = abs(ra - rt), abs(ca - ct)

        # Check for component moves
        if 'R' in piece_moves and (dr == 0 or dc == 0):
            return True
        if 'B' in piece_moves and (dr == dc):
            return True
        if 'N' in piece_moves and ((dr == 1 and dc == 2) or (dr == 2 and dc == 1)):
            return True
        if 'K' in piece_moves and (dr <= 1 and dc <= 1):
            return True
        return False

    # Helper function to find the canonical representation of a position.
    def get_canonical(pos_tuple):
        ra, ca, rk, ck = pos_tuple
        symmetries = set()
        
        # All 8 symmetries of a square board
        transforms = [
            lambda r, c: (r, c),       # Identity
            lambda r, c: (c, 7-r),     # Rotate 90
            lambda r, c: (7-r, 7-c),   # Rotate 180
            lambda r, c: (7-c, r),     # Rotate 270
            lambda r, c: (r, 7-c),     # Flip Horizontal
            lambda r, c: (7-r, c),     # Flip Vertical
            lambda r, c: (c, r),       # Flip Main Diagonal
            lambda r, c: (7-c, 7-r),   # Flip Anti-Diagonal
        ]

        for t in transforms:
            p_a = t(ra, ca)
            p_k = t(rk, ck)
            symmetries.add((p_a[0], p_a[1], p_k[0], p_k[1]))
            
        return min(symmetries)

    # The 7 unique super-pieces formed by combining two distinct standard pieces.
    PIECES = collections.OrderedDict([
        ("Queen (R+B)", frozenset(['R', 'B'])),
        ("Empress (R+N)", frozenset(['R', 'N'])),
        ("Archbishop (B+N)", frozenset(['B', 'N'])),
        ("Amazon (Q+N)", frozenset(['R', 'B', 'N'])),
        ("Dragon King (K+R)", frozenset(['K', 'R'])),
        ("Dragon Horse (K+B)", frozenset(['K', 'B'])),
        ("Knight-King (K+N)", frozenset(['K', 'N']))
    ])

    piece_mate_counts = []
    board_range = range(8)

    for name, moves in PIECES.items():
        unique_mates = set()
        for r_a in board_range:
            for c_a in board_range:
                pos_a = (r_a, c_a)
                for r_k in board_range:
                    for c_k in board_range:
                        if r_a == r_k and c_a == c_k:
                            continue
                        pos_k = (r_k, c_k)

                        # Condition 1: King must be in check.
                        if not is_attacked(pos_a, pos_k, moves):
                            continue

                        # Condition 2: King must have no legal escape moves.
                        has_escape = False
                        for dr_k in [-1, 0, 1]:
                            for dc_k in [-1, 0, 1]:
                                if dr_k == 0 and dc_k == 0:
                                    continue
                                
                                r_e, c_e = r_k + dr_k, c_k + dc_k
                                pos_escape = (r_e, c_e)

                                # Escape must be on the board
                                if not (0 <= r_e < 8 and 0 <= c_e < 8):
                                    continue
                                
                                # Escape must not be attacked or on the attacker's square
                                if pos_escape != pos_a and not is_attacked(pos_a, pos_escape, moves):
                                    has_escape = True
                                    break
                            if has_escape:
                                break
                        
                        if not has_escape:
                            # This is a checkmate position.
                            canonical_mate = get_canonical((r_a, c_a, r_k, c_k))
                            unique_mates.add(canonical_mate)
        
        piece_mate_counts.append(len(unique_mates))
    
    total_checkmates = sum(piece_mate_counts)
    equation_str = " + ".join(map(str, piece_mate_counts)) + f" = {total_checkmates}"
    print(equation_str)

solve_chess_puzzle()
<<<940>>>