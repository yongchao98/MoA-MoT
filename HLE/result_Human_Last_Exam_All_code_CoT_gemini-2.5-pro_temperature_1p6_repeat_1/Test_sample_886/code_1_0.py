import collections

def is_on_board(r, c):
    """Checks if a given coordinate (r, c) is on the 8x8 board."""
    return 0 <= r < 8 and 0 <= c < 8

def get_king_moves(r, c):
    """Returns a set of all valid squares a king can move to from (r, c)."""
    moves = set()
    for dr in [-1, 0, 1]:
        for dc in [-1, 0, 1]:
            if dr == 0 and dc == 0:
                continue
            nr, nc = r + dr, c + dc
            if is_on_board(nr, nc):
                moves.add((nr, nc))
    return moves

# --- Piece Attack Functions ---
# These functions check if a piece at (fr, fc) attacks a target square (tr, tc)
# on an empty board.

def is_pawn_attack(tr, tc, fr, fc):
    """A standard white pawn's forward diagonal attack."""
    return tr - fr == 1 and abs(tc - fc) == 1

def is_knight_attack(tr, tc, fr, fc):
    dr = abs(tr - fr)
    dc = abs(tc - fc)
    return (dr == 1 and dc == 2) or (dr == 2 and dc == 1)

def is_bishop_attack(tr, tc, fr, fc):
    return abs(tr - fr) == abs(tc - fc)

def is_rook_attack(tr, tc, fr, fc):
    return tr == fr or tc == fc

def is_queen_attack(tr, tc, fr, fc):
    return is_rook_attack(tr, tc, fr, fc) or is_bishop_attack(tr, tc, fr, fc)

def is_king_attack(tr, tc, fr, fc):
    return abs(tr - fr) <= 1 and abs(tc - fc) <= 1

def solve():
    """
    Calculates the total number of distinct checkmate positions achievable
    by a hybrid piece (combination of two standard pieces) against a lone king.
    """
    attack_functions = {
        'P': is_pawn_attack, 'N': is_knight_attack, 'B': is_bishop_attack,
        'R': is_rook_attack, 'Q': is_queen_attack, 'K': is_king_attack,
    }
    pieces = list(attack_functions.keys())
    
    piece_pairs = []
    for i in range(len(pieces)):
        for j in range(i + 1, len(pieces)):
            piece_pairs.append((pieces[i], pieces[j]))

    checkmate_positions = set()

    for kr in range(8):
        for kc in range(8):
            k_pos = (kr, kc)
            squares_to_control = get_king_moves(kr, kc)
            squares_to_control.add(k_pos)

            for wr in range(8):
                for wc in range(8):
                    w_pos = (wr, wc)

                    if k_pos == w_pos:
                        continue
                    
                    # Condition: King cannot capture the attacking piece
                    if is_king_attack(wr, wc, kr, kc):
                        continue

                    # Pre-calculate which squares each standard piece attacks
                    coverage = collections.defaultdict(set)
                    for tr, tc in squares_to_control:
                        for piece, attack_func in attack_functions.items():
                            if attack_func(tr, tc, wr, wc):
                                coverage[piece].add((tr, tc))
                    
                    # Check if any hybrid piece delivers checkmate
                    for p1, p2 in piece_pairs:
                        combined_coverage = coverage[p1].union(coverage[p2])
                        
                        is_in_check = k_pos in combined_coverage
                        all_escapes_covered = squares_to_control.issubset(combined_coverage)

                        if is_in_check and all_escapes_covered:
                            checkmate_positions.add((w_pos, k_pos))
                            # Position is a mate, no need to check other pairs
                            break 
                            
    final_count = len(checkmate_positions)
    print("The final number of distinct checkmate positions is:")
    # Per the instructions, we output the number involved in the final result.
    print(final_count)

solve()