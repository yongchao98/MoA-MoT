import itertools

def solve_chess_mates():
    """
    Calculates the number of legal checkmate positions with White (King, Knight, Pawn)
    vs Black (King).
    """

    # Helper functions to convert between square index (0-63) and (row, col) coordinates
    def get_coords(sq):
        return sq // 8, sq % 8

    def get_sq(row, col):
        return row * 8 + col

    # Pre-compute attack sets for all pieces on all squares to speed up the main loop
    KING_ATTACKS = {sq: set() for sq in range(64)}
    KNIGHT_ATTACKS = {sq: set() for sq in range(64)}
    PAWN_ATTACKS = {sq: set() for sq in range(64)}

    for r in range(8):
        for c in range(8):
            sq = get_sq(r, c)
            
            # King attacks
            for dr in [-1, 0, 1]:
                for dc in [-1, 0, 1]:
                    if dr == 0 and dc == 0: continue
                    nr, nc = r + dr, c + dc
                    if 0 <= nr < 8 and 0 <= nc < 8:
                        KING_ATTACKS[sq].add(get_sq(nr, nc))

            # Knight attacks
            for dr, dc in [(-2, -1), (-2, 1), (-1, -2), (-1, 2), (1, -2), (1, 2), (2, -1), (2, 1)]:
                nr, nc = r + dr, c + dc
                if 0 <= nr < 8 and 0 <= nc < 8:
                    KNIGHT_ATTACKS[sq].add(get_sq(nr, nc))
            
            # White Pawn attacks
            nr = r + 1 # Pawns move to higher ranks
            if 0 <= nr < 8:
                if c > 0: PAWN_ATTACKS[sq].add(get_sq(nr, c - 1))
                if c < 7: PAWN_ATTACKS[sq].add(get_sq(nr, c + 1))

    # Initialize counters
    knight_mates = 0
    pawn_mates = 0
    double_check_mates = 0

    # Iterate through all permutations of 4 distinct squares for the 4 pieces
    # wk: White King, bk: Black King, wn: White Knight, wp: White Pawn
    for wk_sq, bk_sq, wn_sq, wp_sq in itertools.permutations(range(64), 4):
        
        # 1. Legality Check: King Proximity
        # The two kings cannot be on adjacent squares (Chebyshev distance > 1)
        r1, c1 = get_coords(wk_sq)
        r2, c2 = get_coords(bk_sq)
        if max(abs(r1 - r2), abs(c1 - c2)) <= 1:
            continue

        # 2. Legality Check: Pawn Placement
        # White pawns cannot be on the 1st (rank 0) or 8th (rank 7)
        if get_coords(wp_sq)[0] in [0, 7]:
            continue

        # 3. Checkmate Condition: Is Black in check?
        knight_attack_set = KNIGHT_ATTACKS[wn_sq]
        pawn_attack_set = PAWN_ATTACKS[wp_sq]
        
        checked_by_knight = bk_sq in knight_attack_set
        checked_by_pawn = bk_sq in pawn_attack_set

        if not checked_by_knight and not checked_by_pawn:
            continue
            
        # 4. Checkmate Condition: Can the Black King escape?
        # A move is legal if the destination square is not attacked by any white piece.
        
        # Get all squares controlled by White
        all_white_attacks = KING_ATTACKS[wk_sq].union(knight_attack_set).union(pawn_attack_set)
        
        # Check every potential escape square for the black king
        has_legal_move = False
        for move_sq in KING_ATTACKS[bk_sq]:
            if move_sq not in all_white_attacks:
                has_legal_move = True
                break
        
        # If the king is in check and has no legal moves, it's checkmate
        if not has_legal_move:
            if checked_by_knight and checked_by_pawn:
                double_check_mates += 1
            elif checked_by_knight:
                knight_mates += 1
            else: # checked_by_pawn
                pawn_mates += 1

    total_count = knight_mates + pawn_mates + double_check_mates
    
    # Print the results in the requested equation format
    print(f"{knight_mates} + {pawn_mates} + {double_check_mates} = {total_count}")
    print(f"In total, there are {total_count} possible checkmate positions.")

solve_chess_mates()