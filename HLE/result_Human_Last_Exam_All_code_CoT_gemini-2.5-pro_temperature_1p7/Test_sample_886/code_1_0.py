import functools

def solve():
    """
    This script calculates the total number of distinct checkmate positions achievable 
    with a single hybrid chess piece against a lone king. A hybrid piece is formed 
    by combining the movement capabilities of two distinct standard chess pieces 
    (pawn excluded due to ambiguity).
    """

    # --- Step 1: Define base piece movements with caching for performance ---
    @functools.lru_cache(maxsize=None)
    def get_king_attacks(pos):
        r, c = pos
        attacks = set()
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0:
                    continue
                nr, nc = r + dr, c + dc
                if 0 <= nr < 8 and 0 <= nc < 8:
                    attacks.add((nr, nc))
        return attacks

    @functools.lru_cache(maxsize=None)
    def get_knight_attacks(pos):
        r, c = pos
        attacks = set()
        moves = [(-2, -1), (-2, 1), (-1, -2), (-1, 2),
                 (1, -2), (1, 2), (2, -1), (2, 1)]
        for dr, dc in moves:
            nr, nc = r + dr, c + dc
            if 0 <= nr < 8 and 0 <= nc < 8:
                attacks.add((nr, nc))
        return attacks

    @functools.lru_cache(maxsize=None)
    def get_bishop_attacks(pos):
        r, c = pos
        attacks = set()
        for dr, dc in [(-1, -1), (-1, 1), (1, -1), (1, 1)]:
            for i in range(1, 8):
                nr, nc = r + i * dr, c + i * dc
                if 0 <= nr < 8 and 0 <= nc < 8:
                    attacks.add((nr, nc))
                else:
                    break
        return attacks

    @functools.lru_cache(maxsize=None)
    def get_rook_attacks(pos):
        r, c = pos
        attacks = set()
        for i in range(8):
            if i != r:
                attacks.add((i, c))
            if i != c:
                attacks.add((r, i))
        return attacks

    # --- Step 2: Define hybrid pieces by combining move sets ---
    HYBRID_PIECES = {
        "Amazon (Queen + Knight)": (get_rook_attacks, get_bishop_attacks, get_knight_attacks),
        "Chancellor (Rook + Knight)": (get_rook_attacks, get_knight_attacks),
        "Archbishop (Bishop + Knight)": (get_bishop_attacks, get_knight_attacks),
        "King + Knight": (get_king_attacks, get_knight_attacks),
        "King + Rook": (get_king_attacks, get_rook_attacks),
        "King + Bishop": (get_king_attacks, get_bishop_attacks),
    }

    @functools.lru_cache(maxsize=None)
    def get_hybrid_attacks(pos, move_funcs_tuple):
        """Computes all attacks for a hybrid piece by combining base move sets."""
        all_attacks = set()
        for move_func in move_funcs_tuple:
            all_attacks.update(move_func(pos))
        return all_attacks

    # --- Step 3: Implement the checkmate logic ---
    def is_checkmate(king_pos, piece_pos, move_funcs_tuple):
        """Checks if a king is checkmated by a single hybrid piece."""
        # Condition 3: King can't capture the piece (it's not adjacent).
        if piece_pos in get_king_attacks(king_pos):
            return False

        piece_attacks = get_hybrid_attacks(piece_pos, move_funcs_tuple)

        # Condition 1: King must be in check.
        if king_pos not in piece_attacks:
            return False

        # Condition 2: All king's escape squares must be attacked.
        for escape_square in get_king_attacks(king_pos):
            if escape_square not in piece_attacks:
                return False

        return True

    # --- Step 4 & 5: Count all positions and display results ---
    all_squares = [(r, c) for r in range(8) for c in range(8)]
    total_checkmates = 0
    piece_mate_counts = []
    
    print("Calculating distinct checkmate positions for each hybrid piece...")
    
    for name, move_funcs in HYBRID_PIECES.items():
        mate_count = 0
        for king_pos in all_squares:
            for piece_pos in all_squares:
                if king_pos == piece_pos:
                    continue
                # The tuple of functions must be passed to the cached is_checkmate
                if is_checkmate(king_pos, piece_pos, move_funcs):
                    mate_count += 1
        
        piece_mate_counts.append((name, mate_count))
        total_checkmates += mate_count
    
    print("\n--- Results ---")
    equation_parts = []
    for name, count in piece_mate_counts:
        print(f"Checkmates for {name}: {count}")
        equation_parts.append(str(count))
        
    print("\nTotal distinct checkmate positions:")
    equation = " + ".join(equation_parts)
    print(f"{equation} = {total_checkmates}")
    
    print(f"\n<<<{total_checkmates}>>>")

solve()