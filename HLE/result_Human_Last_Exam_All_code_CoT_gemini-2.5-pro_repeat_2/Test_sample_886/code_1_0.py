import collections

def solve_chess_mate_problem():
    """
    This script calculates the total number of distinct checkmate positions
    achievable on a chessboard with a single hybrid piece against a lone king.
    A hybrid piece is formed by combining the moves of two standard chess pieces.
    """
    BOARD_SIZE = 8

    # Memoization caches for performance
    _attacks_cache = {}

    def is_valid(r, c):
        return 0 <= r < BOARD_SIZE and 0 <= c < BOARD_SIZE

    def get_attacks_from_cache(piece_type, r, c):
        key = (piece_type, r, c)
        if key in _attacks_cache:
            return _attacks_cache[key]
        return None

    def set_attacks_in_cache(piece_type, r, c, attacks):
        key = (piece_type, r, c)
        _attacks_cache[key] = attacks

    def get_rook_attacks(r, c):
        cached = get_attacks_from_cache('R', r, c)
        if cached is not None: return cached
        attacks = set()
        for i in range(BOARD_SIZE):
            if i != r: attacks.add((i, c))
            if i != c: attacks.add((r, i))
        set_attacks_in_cache('R', r, c, attacks)
        return attacks

    def get_bishop_attacks(r, c):
        cached = get_attacks_from_cache('B', r, c)
        if cached is not None: return cached
        attacks = set()
        for dr, dc in [(1, 1), (1, -1), (-1, 1), (-1, -1)]:
            for i in range(1, BOARD_SIZE):
                nr, nc = r + i * dr, c + i * dc
                if is_valid(nr, nc):
                    attacks.add((nr, nc))
                else:
                    break
        set_attacks_in_cache('B', r, c, attacks)
        return attacks

    def get_knight_attacks(r, c):
        cached = get_attacks_from_cache('N', r, c)
        if cached is not None: return cached
        attacks = set()
        moves = [(1, 2), (1, -2), (-1, 2), (-1, -2), (2, 1), (2, -1), (-2, 1), (-2, -1)]
        for dr, dc in moves:
            nr, nc = r + dr, c + dc
            if is_valid(nr, nc):
                attacks.add((nr, nc))
        set_attacks_in_cache('N', r, c, attacks)
        return attacks

    def get_king_attacks(r, c):
        cached = get_attacks_from_cache('K', r, c)
        if cached is not None: return cached
        attacks = set()
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0: continue
                nr, nc = r + dr, c + dc
                if is_valid(nr, nc): attacks.add((nr, nc))
        set_attacks_in_cache('K', r, c, attacks)
        return attacks

    def get_pawn_attacks(r, c):
        cached = get_attacks_from_cache('P', r, c)
        if cached is not None: return cached
        # An omni-directional pawn attacks the four adjacent diagonal squares.
        attacks = set()
        for dr, dc in [(1, 1), (1, -1), (-1, 1), (-1, -1)]:
            nr, nc = r + dr, c + dc
            if is_valid(nr, nc): attacks.add((nr, nc))
        set_attacks_in_cache('P', r, c, attacks)
        return attacks

    def get_queen_attacks(r, c):
        return get_rook_attacks(r, c).union(get_bishop_attacks(r, c))

    PIECE_ATTACKS = {'R': get_rook_attacks, 'B': get_bishop_attacks, 'N': get_knight_attacks,
                     'Q': get_queen_attacks, 'K': get_king_attacks, 'P': get_pawn_attacks}

    def get_hybrid_attacks(piece_types, r, c):
        attacks = set()
        for p_type in piece_types:
            attacks.update(PIECE_ATTACKS[p_type](r, c))
        return attacks

    def check_mate(king_pos, hybrid_pos, hybrid_piece_types):
        king_escape_squares = get_king_attacks(king_pos[0], king_pos[1])
        if hybrid_pos in king_escape_squares:
            return False
        
        hybrid_attack_squares = get_hybrid_attacks(hybrid_piece_types, hybrid_pos[0], hybrid_pos[1])
        if king_pos not in hybrid_attack_squares:
            return False
        
        return king_escape_squares.issubset(hybrid_attack_squares)

    def get_canonical(k_pos, h_pos):
        positions = []
        r_k, c_k = k_pos
        r_h, c_h = h_pos
        for _ in range(4):
            # Original and horizontal flip
            positions.append(tuple(sorted(((r_k, c_k), (r_h, c_h)))))
            positions.append(tuple(sorted(((r_k, 7 - c_k), (r_h, 7 - c_h)))))
            # Rotate 90 degrees: (r, c) -> (c, 7-r)
            r_k, c_k = c_k, 7 - r_k
            r_h, c_h = c_h, 7 - r_h
        return min(positions)

    pieces = ['P', 'N', 'B', 'R', 'Q', 'K']
    combinations = []
    for i in range(len(pieces)):
        for j in range(i + 1, len(pieces)):
            combinations.append((pieces[i], pieces[j]))

    unique_hybrid_movesets = collections.defaultdict(list)
    for p1, p2 in combinations:
        p1_comp = {'R', 'B'} if p1 == 'Q' else {p1}
        p2_comp = {'R', 'B'} if p2 == 'Q' else {p2}
        components = tuple(sorted(list(p1_comp.union(p2_comp))))
        unique_hybrid_movesets[components].append(f"{p1}+{p2}")
    
    sorted_movesets = sorted(list(unique_hybrid_movesets.keys()))
    
    all_found_mates = set()
    final_equation_parts = []
    contributing_pieces = []
    
    for components in sorted_movesets:
        mates_for_this_piece = set()
        # This check can be skipped for performance, but we run it to be exhaustive
        # if 'K' in components: continue 
        
        for r_k in range(BOARD_SIZE):
            for c_k in range(BOARD_SIZE):
                for r_h in range(BOARD_SIZE):
                    for c_h in range(BOARD_SIZE):
                        if (r_k, c_k) == (r_h, c_h): continue
                        if check_mate((r_k, c_k), (r_h, c_h), components):
                            mates_for_this_piece.add(get_canonical((r_k, c_k), (r_h, c_h)))
        
        if len(mates_for_this_piece) > 0:
            names = unique_hybrid_movesets[components]
            rep_name = min(names, key=len)
            contributing_pieces.append(rep_name)
            final_equation_parts.append(str(len(mates_for_this_piece)))

        all_found_mates.update(mates_for_this_piece)

    total_distinct_positions = len(all_found_mates)

    print("The number of distinct checkmate positions for each successful hybrid piece combination are:")
    for piece, count in zip(contributing_pieces, final_equation_parts):
        print(f"{count}", end='')
        print(f" positions from the {piece} piece")

    print("\nAfter accounting for overlaps (where the same checkmate position can be achieved by different pieces), the final result is calculated.")
    print(f"Total distinct checkmate positions: {total_distinct_positions}")

solve_chess_mate_problem()