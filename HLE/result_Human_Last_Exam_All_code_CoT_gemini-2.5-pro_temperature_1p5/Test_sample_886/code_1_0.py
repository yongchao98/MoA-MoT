def solve_chess_puzzle():
    """
    This script calculates the total number of distinct checkmate positions achievable
    by a single hybrid chess piece against a lone king.
    A hybrid piece is one that combines the move sets of two distinct standard pieces.
    """

    # Global dictionaries to store precomputed moves for performance
    ROOK_ATTACKS = [set() for _ in range(64)]
    BISHOP_ATTACKS = [set() for _ in range(64)]
    KNIGHT_ATTACKS = [set() for _ in range(64)]
    KING_MOVES = [set() for _ in range(64)]
    PAWN_ATTACKS = [set() for _ in range(64)]  # Symmetrical pawn attacks

    def is_valid(r, c):
        """Check if a square (r, c) is on the 8x8 board."""
        return 0 <= r < 8 and 0 <= c < 8

    def precompute_all_moves():
        """Fills the global move dictionaries for each square on the board."""
        for r_start in range(8):
            for c_start in range(8):
                pos = r_start * 8 + c_start

                # Rook moves
                for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                    r, c = r_start + dr, c_start + dc
                    while is_valid(r, c):
                        ROOK_ATTACKS[pos].add(r * 8 + c)
                        r, c = r + dr, c + dc

                # Bishop moves
                for dr, dc in [(1, 1), (1, -1), (-1, 1), (-1, -1)]:
                    r, c = r_start + dr, c_start + dc
                    while is_valid(r, c):
                        BISHOP_ATTACKS[pos].add(r * 8 + c)
                        r, c = r + dr, c + dc

                # Knight moves
                for dr, dc in [(1, 2), (1, -2), (-1, 2), (-1, -2),
                               (2, 1), (2, -1), (-2, 1), (-2, -1)]:
                    r, c = r_start + dr, c_start + dc
                    if is_valid(r, c):
                        KNIGHT_ATTACKS[pos].add(r * 8 + c)

                # King moves
                for dr in [-1, 0, 1]:
                    for dc in [-1, 0, 1]:
                        if dr == 0 and dc == 0:
                            continue
                        r, c = r_start + dr, c_start + dc
                        if is_valid(r, c):
                            KING_MOVES[pos].add(r * 8 + c)

                # Symmetrical Pawn attacks (one square diagonally in 4 directions)
                for dr, dc in [(1, 1), (1, -1), (-1, 1), (-1, -1)]:
                    r, c = r_start + dr, c_start + dc
                    if is_valid(r, c):
                        PAWN_ATTACKS[pos].add(r * 8 + c)

    def get_hybrid_moves(h_pos, move_sets):
        """Get all squares a hybrid piece can move to/attack."""
        moves = set()
        if 'R' in move_sets: moves.update(ROOK_ATTACKS[h_pos])
        if 'B' in move_sets: moves.update(BISHOP_ATTACKS[h_pos])
        if 'N' in move_sets: moves.update(KNIGHT_ATTACKS[h_pos])
        if 'K' in move_sets: moves.update(KING_MOVES[h_pos])
        if 'P' in move_sets: moves.update(PAWN_ATTACKS[h_pos])
        return moves

    def is_checkmate(h_pos, k_pos, move_sets):
        """
        Determines if a king at k_pos is checkmated by a hybrid piece at h_pos.
        """
        has_king_move = 'K' in move_sets
        
        # An illegal position: a king cannot be adjacent to the opposing king.
        # If the hybrid has king moves, it projects this "royal zone".
        # Also, check cannot be delivered by a king's move.
        if has_king_move and k_pos in KING_MOVES[h_pos]:
            return False

        # Determine which moves can deliver check (non-king moves)
        check_delivering_moves = get_hybrid_moves(h_pos, move_sets - {'K'})

        # Condition 1: The king must be in check from a non-king move.
        if k_pos not in check_delivering_moves:
            return False

        all_hybrid_attacks = get_hybrid_moves(h_pos, move_sets)

        # Condition 2: All king's escape squares must be controlled.
        for escape_pos in KING_MOVES[k_pos]:
            is_attacked = escape_pos in all_hybrid_attacks
            is_adjacent_to_hybrid_king = has_king_move and (escape_pos in KING_MOVES[h_pos])

            if not is_attacked and not is_adjacent_to_hybrid_king:
                return False  # Found a legal move, so not checkmate.

        return True  # King is in check and has no legal moves.

    precompute_all_moves()

    # Define the movement combinations for the hybrid pieces that can checkmate.
    mating_hybrids = {
        "Queen (Rook+Bishop)": frozenset({'R', 'B'}),
        "Chancellor (Rook+Knight)": frozenset({'R', 'N'}),
        "Archbishop (Bishop+Knight)": frozenset({'B', 'N'}),
        "Rook+King": frozenset({'R', 'K'}),
        "Knight+King": frozenset({'N', 'K'}),
        "Rook+Pawn": frozenset({'R', 'P'})
    }

    total_checkmate_positions = 0
    equation_parts = []

    for name, move_sets in mating_hybrids.items():
        mate_positions_for_piece = set()
        for k_pos in range(64):
            for h_pos in range(64):
                if k_pos == h_pos:
                    continue

                if is_checkmate(h_pos, k_pos, move_sets):
                    mate_positions_for_piece.add((k_pos, h_pos))
        
        num_mates = len(mate_positions_for_piece)
        equation_parts.append(str(num_mates))
        total_checkmate_positions += num_mates

    final_equation = " + ".join(equation_parts) + f" = {total_checkmate_positions}"
    print(final_equation)
    print(f"<<<{total_checkmate_positions}>>>")

solve_chess_puzzle()