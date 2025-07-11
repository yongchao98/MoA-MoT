import collections

def solve():
    """
    Calculates the number of distinct checkmate positions on a chessboard 
    with a single super-piece against a lone king.
    """

    # 1. Define movement generation functions for standard pieces
    def get_rook_moves(r, c):
        moves = set()
        for i in range(8):
            if i != r: moves.add((i, c))
            if i != c: moves.add((r, i))
        return moves

    def get_bishop_moves(r, c):
        moves = set()
        for i in range(1, 8):
            for dr, dc in [(1, 1), (1, -1), (-1, 1), (-1, -1)]:
                nr, nc = r + i * dr, c + i * dc
                if 0 <= nr < 8 and 0 <= nc < 8:
                    moves.add((nr, nc))
        return moves

    def get_knight_moves(r, c):
        moves = set()
        for dr, dc in [(-2, -1), (-2, 1), (-1, -2), (-1, 2), 
                       (1, -2), (1, 2), (2, -1), (2, 1)]:
            nr, nc = r + dr, c + dc
            if 0 <= nr < 8 and 0 <= nc < 8:
                moves.add((nr, nc))
        return moves

    def get_king_moves(r, c):
        moves = set()
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0: continue
                nr, nc = r + dr, c + dc
                if 0 <= nr < 8 and 0 <= nc < 8:
                    moves.add((nr, nc))
        return moves

    def get_queen_moves(r, c):
        return get_rook_moves(r, c).union(get_bishop_moves(r, c))

    # 2. Define the 7 unique super-pieces from combinations
    PIECE_MOVE_GENERATORS = collections.OrderedDict([
        ("Queen (Rook+Bishop)", [get_rook_moves, get_bishop_moves]),
        ("Chancellor (Rook+Knight)", [get_rook_moves, get_knight_moves]),
        ("Gryphon (Rook+King)", [get_rook_moves, get_king_moves]),
        ("Archbishop (Bishop+Knight)", [get_bishop_moves, get_knight_moves]),
        ("Centaur (Bishop+King)", [get_bishop_moves, get_king_moves]),
        ("Amazon (Queen+Knight)", [get_queen_moves, get_knight_moves]),
        ("Kight (Knight+King)", [get_knight_moves, get_king_moves]),
    ])

    # Helper function to get all attacked squares for a piece
    def get_attack_set(piece_pos, move_funcs):
        r, c = piece_pos
        attacks = set()
        for func in move_funcs:
            attacks.update(func(r, c))
        return attacks

    # Helper function to normalize a position using board symmetries
    def normalize_position(s_pos, k_pos):
        s_r, s_c = s_pos
        k_r, k_c = k_pos
        
        positions = []
        # Generate 8 symmetries (rotations and reflections)
        for _ in range(2): # Original and horizontal reflection
            current_sr, current_sc = s_r, s_c
            current_kr, current_kc = k_r, k_c
            for _ in range(4): # 4 rotations
                positions.append(tuple(sorted(((current_sr, current_sc), (current_kr, current_kc)))))
                # Rotate 90 degrees: (r, c) -> (c, 7-r)
                current_sr, current_sc = current_sc, 7 - current_sr
                current_kr, current_kc = current_kc, 7 - current_kr
            # Flip horizontally for next iteration: (r, c) -> (r, 7-c)
            s_c, k_c = 7 - s_c, 7 - k_c

        return min(positions)
    
    # 3. Main loop to find and count all mate positions
    master_mate_set = set()
    piece_mate_counts = {}
    
    print("Calculating distinct checkmate positions for each piece type...\n")
    
    for piece_name, move_funcs in PIECE_MOVE_GENERATORS.items():
        current_piece_mates = set()
        for sr in range(8):
            for sc in range(8):
                s_pos = (sr, sc)
                attack_set = get_attack_set(s_pos, move_funcs)
                
                for kr in range(8):
                    for kc in range(8):
                        k_pos = (kr, kc)
                        
                        if s_pos == k_pos: continue

                        # Rule 1: King must be in check
                        if k_pos not in attack_set:
                            continue
                            
                        # Rule 2: Super-piece cannot be adjacent to King
                        if max(abs(sr - kr), abs(sc - kc)) == 1:
                            continue
                            
                        # Rule 3: All King's escape squares must be attacked
                        is_mate = True
                        escape_squares = get_king_moves(kr, kc)
                        for esc_pos in escape_squares:
                            if esc_pos not in attack_set:
                                is_mate = False
                                break
                        
                        if is_mate:
                            canonical_form = normalize_position(s_pos, k_pos)
                            current_piece_mates.add(canonical_form)
        
        piece_mate_counts[piece_name] = len(current_piece_mates)
        master_mate_set.update(current_piece_mates)

    # 4. Output the results
    equation_parts = []
    for piece_name, count in piece_mate_counts.items():
        print(f"Number of distinct checkmate positions for {piece_name}: {count}")
        equation_parts.append(str(count))
    
    final_total = len(master_mate_set)
    
    # Building a representative final equation string
    equation_str = " + ".join(equation_parts)

    print(f"\nThe individual counts are: {equation_str}")
    print(f"(Note: The sum of individual counts is not the total, as positions can be created by multiple pieces)")
    print("\n---")
    print("The final total number of distinct checkmate positions from any combined piece is:")
    print(final_total)

solve()
<<<944>>>