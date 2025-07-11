import collections

def solve():
    """
    Calculates the total number of distinct checkmate positions achievable
    by a single piece that combines the moves of two standard chess pieces.
    """

    # Helper functions to convert between (row, col) and square index
    def to_sq(r, c):
        if 0 <= r < 8 and 0 <= c < 8:
            return r * 8 + c
        return -1

    # --- Pre-calculate attack patterns for each piece from each square ---

    KING_ATTACKS = []
    for r in range(8):
        for c in range(8):
            attacks = set()
            for dr in [-1, 0, 1]:
                for dc in [-1, 0, 1]:
                    if dr == 0 and dc == 0:
                        continue
                    sq = to_sq(r + dr, c + dc)
                    if sq != -1:
                        attacks.add(sq)
            KING_ATTACKS.append(attacks)

    KNIGHT_ATTACKS = []
    for r in range(8):
        for c in range(8):
            attacks = set()
            for dr, dc in [(1, 2), (1, -2), (-1, 2), (-1, -2), (2, 1), (2, -1), (-2, 1), (-2, -1)]:
                sq = to_sq(r + dr, c + dc)
                if sq != -1:
                    attacks.add(sq)
            KNIGHT_ATTACKS.append(attacks)

    BISHOP_ATTACKS = []
    for r_start in range(8):
        for c_start in range(8):
            attacks = set()
            for dr, dc in [(-1, -1), (-1, 1), (1, -1), (1, 1)]:
                for i in range(1, 8):
                    r, c = r_start + i * dr, c_start + i * dc
                    if 0 <= r < 8 and 0 <= c < 8:
                        attacks.add(to_sq(r, c))
                    else:
                        break
            BISHOP_ATTACKS.append(attacks)

    ROOK_ATTACKS = []
    for r_start in range(8):
        for c_start in range(8):
            attacks = set()
            for dr, dc in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
                for i in range(1, 8):
                    r, c = r_start + i * dr, c_start + i * dc
                    if 0 <= r < 8 and 0 <= c < 8:
                        attacks.add(to_sq(r, c))
                    else:
                        break
            ROOK_ATTACKS.append(attacks)

    QUEEN_ATTACKS = [b | r for b, r in zip(BISHOP_ATTACKS, ROOK_ATTACKS)]

    # Pawn attacks are its diagonal captures. Assuming a 'white' piece moving up the board.
    PAWN_ATTACKS = []
    for r in range(8):
        for c in range(8):
            attacks = set()
            if r > 0: # Can't attack from the top rank
                if c > 0: attacks.add(to_sq(r - 1, c - 1))
                if c < 7: attacks.add(to_sq(r - 1, c + 1))
            PAWN_ATTACKS.append(attacks)

    # --- Define the unique super-pieces ---

    PIECES = ["P", "N", "B", "R", "Q", "K"]
    ATTACKS_MAP = {
        "P": PAWN_ATTACKS, "N": KNIGHT_ATTACKS, "B": BISHOP_ATTACKS,
        "R": ROOK_ATTACKS, "Q": QUEEN_ATTACKS, "K": KING_ATTACKS
    }

    # Generate unique super-pieces based on their combined move sets
    super_pieces = collections.OrderedDict()
    for i in range(len(PIECES)):
        for j in range(i + 1, len(PIECES)):
            p1_name, p2_name = PIECES[i], PIECES[j]
            
            # Simplify combinations involving Queen (Q = B+R) to get unique move sets
            moves = {p1_name, p2_name}
            if 'Q' in moves:
                if 'B' in moves: moves.remove('B')
                if 'R' in moves: moves.remove('R')
            if 'B' in moves and 'R' in moves:
                moves = {'Q'}
            
            name = "+".join(sorted(list(moves)))
            if name not in super_pieces:
                super_pieces[name] = [ATTACKS_MAP[p] for p in sorted(list(moves))]

    # --- Calculate checkmate positions for each super-piece ---

    piece_mate_counts = collections.OrderedDict()

    for name, attack_sets in super_pieces.items():
        count = 0
        is_king_combo = any(s == KING_ATTACKS for s in attack_sets)

        for wp_sq in range(64): # White super-piece square
            for bk_sq in range(64): # Black king square
                if wp_sq == bk_sq:
                    continue

                # 1. Get all squares attacked by the super-piece
                wp_attacked_squares = set()
                for attack_set in attack_sets:
                    wp_attacked_squares.update(attack_set[wp_sq])

                # 2. King must be in check
                if bk_sq not in wp_attacked_squares:
                    continue

                # 3. King must have no legal moves
                has_legal_move = False
                bk_potential_moves = KING_ATTACKS[bk_sq]
                
                # Check escape to empty squares
                for move_sq in bk_potential_moves:
                    if move_sq != wp_sq and move_sq not in wp_attacked_squares:
                        has_legal_move = True
                        break
                if has_legal_move:
                    continue

                # Check capture of the checking piece
                if wp_sq in bk_potential_moves:
                    # King can capture if the piece's square is not self-defended.
                    # Only king-combos can defend their own square.
                    if not is_king_combo:
                        has_legal_move = True
                
                if not has_legal_move:
                    count += 1
        
        piece_mate_counts[name] = count

    # --- Print the results ---
    
    total_mates = 0
    equation_parts = []
    
    # Sort keys for consistent output, with single letters first
    sorted_keys = sorted(piece_mate_counts.keys(), key=lambda x: (len(x), x))

    print("Number of distinct checkmate positions for each unique piece combination:")
    for name in sorted_keys:
        count = piece_mate_counts[name]
        print(f"{name}: {count}")
        total_mates += count
        equation_parts.append(str(count))
        
    print("\nTotal checkmate positions:")
    print(f"{' + '.join(equation_parts)} = {total_mates}")

solve()
<<<856>>>