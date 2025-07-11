import itertools

def solve_chess_checkmates():
    """
    Calculates the number of legal checkmate positions with White (King, Knight, Pawn)
    vs Black (King) by iterating through all possible piece placements.
    
    This is a computationally intensive task and may take a few minutes to run.
    """
    
    # Pre-compute attack maps for all pieces from all squares for efficiency.
    king_attack_map = [set() for _ in range(64)]
    knight_attack_map = [set() for _ in range(64)]
    pawn_attack_map = [set() for _ in range(64)]

    for sq in range(64):
        r, c = sq // 8, sq % 8
        
        # King attacks
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0:
                    continue
                nr, nc = r + dr, c + dc
                if 0 <= nr < 8 and 0 <= nc < 8:
                    king_attack_map[sq].add(nr * 8 + nc)

        # Knight attacks
        for dr, dc in [(-2, -1), (-2, 1), (-1, -2), (-1, 2), (1, -2), (1, 2), (2, -1), (2, 1)]:
            nr, nc = r + dr, c + dc
            if 0 <= nr < 8 and 0 <= nc < 8:
                knight_attack_map[sq].add(nr * 8 + nc)

        # White Pawn attacks (moves from low rank to high rank)
        nr = r + 1
        if 0 <= nr < 8:
            if c > 0: pawn_attack_map[sq].add(nr * 8 + c - 1)
            if c < 7: pawn_attack_map[sq].add(nr * 8 + c + 1)
            
    count_by_knight = 0
    count_by_pawn = 0
    count_by_double = 0
    
    squares = range(64)
    # Generate all permutations of 4 pieces on 64 squares.
    for p in itertools.permutations(squares, 4):
        wk_sq, wn_sq, wp_sq, bk_sq = p

        # --- Filter 1: Legality of Piece Placement ---
        
        # A pawn on the 1st or 8th rank is illegal (would have promoted).
        if wp_sq // 8 in {0, 7}:
            continue
            
        # Kings on adjacent squares are illegal. This also ensures White isn't in check.
        if bk_sq in king_attack_map[wk_sq]:
            continue

        # --- Filter 2: Checkmate Conditions ---
        
        # Condition A: Is the Black King in check?
        knight_attacks = knight_attack_map[wn_sq]
        pawn_attacks = pawn_attack_map[wp_sq]
        
        is_check_by_knight = bk_sq in knight_attacks
        is_check_by_pawn = bk_sq in pawn_attacks

        if not (is_check_by_knight or is_check_by_pawn):
            continue

        # Condition B: Can the Black King escape?
        king_attacks = king_attack_map[wk_sq]
        all_white_attacks = king_attacks | knight_attacks | pawn_attacks
        
        black_king_escapes = king_attack_map[bk_sq]

        has_legal_move = False
        for escape_square in black_king_escapes:
            if escape_square not in all_white_attacks:
                has_legal_move = True
                break
        
        if has_legal_move:
            continue

        # If all checks pass, this is a valid checkmate. Categorize and count it.
        if is_check_by_knight and is_check_by_pawn:
            count_by_double += 1
        elif is_check_by_knight:
            count_by_knight += 1
        else: # is_check_by_pawn must be true
            count_by_pawn += 1
            
    total_checkmates = count_by_knight + count_by_pawn + count_by_double

    print("Calculation complete.")
    print(f"Checkmates delivered by Knight: {count_by_knight}")
    print(f"Checkmates delivered by Pawn: {count_by_pawn}")
    print(f"Checkmates delivered by Double Check: {count_by_double}")
    print("\nFinal Equation:")
    print(f"{count_by_knight} + {count_by_pawn} + {count_by_double} = {total_checkmates}")
    
    # Store final answer for parsing
    print(f"\n<<<{total_checkmates}>>>")

solve_chess_checkmates()