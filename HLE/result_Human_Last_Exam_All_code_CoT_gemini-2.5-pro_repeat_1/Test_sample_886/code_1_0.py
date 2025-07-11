def solve_chess_checkmates():
    """
    Calculates the total number of distinct checkmate positions achievable
    by a single hybrid chess piece (Amazon, Chancellor, Archbishop) against a lone king.
    """

    def is_on_board(x, y):
        """Checks if coordinates are within the 8x8 board."""
        return 0 <= x < 8 and 0 <= y < 8

    # --- Attack Pattern Calculation Functions ---

    memoized_attacks = {}

    def get_rook_attacks(x, y):
        """Returns a set of squares attacked by a rook from (x, y)."""
        if ('R', x, y) in memoized_attacks:
            return memoized_attacks[('R', x, y)]
        
        attacks = set()
        for dx, dy in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            for i in range(1, 8):
                nx, ny = x + i * dx, y + i * dy
                if is_on_board(nx, ny):
                    attacks.add((nx, ny))
                else:
                    break
        memoized_attacks[('R', x, y)] = attacks
        return attacks

    def get_bishop_attacks(x, y):
        """Returns a set of squares attacked by a bishop from (x, y)."""
        if ('B', x, y) in memoized_attacks:
            return memoized_attacks[('B', x, y)]

        attacks = set()
        for dx, dy in [(1, 1), (1, -1), (-1, 1), (-1, -1)]:
            for i in range(1, 8):
                nx, ny = x + i * dx, y + i * dy
                if is_on_board(nx, ny):
                    attacks.add((nx, ny))
                else:
                    break
        memoized_attacks[('B', x, y)] = attacks
        return attacks

    def get_knight_attacks(x, y):
        """Returns a set of squares attacked by a knight from (x, y)."""
        if ('N', x, y) in memoized_attacks:
            return memoized_attacks[('N', x, y)]

        attacks = set()
        moves = [(1, 2), (1, -2), (-1, 2), (-1, -2),
                 (2, 1), (2, -1), (-2, 1), (-2, -1)]
        for dx, dy in moves:
            nx, ny = x + dx, y + dy
            if is_on_board(nx, ny):
                attacks.add((nx, ny))
        memoized_attacks[('N', x, y)] = attacks
        return attacks

    # --- Hybrid Piece Definitions ---
    
    hybrid_pieces = {
        "Amazon": lambda x, y: get_rook_attacks(x, y) | get_bishop_attacks(x, y) | get_knight_attacks(x, y),
        "Chancellor": lambda x, y: get_rook_attacks(x, y) | get_knight_attacks(x, y),
        "Archbishop": lambda x, y: get_bishop_attacks(x, y) | get_knight_attacks(x, y),
    }

    mate_counts = {}

    # --- Main Calculation Loop ---
    
    for piece_name, get_piece_attacks in hybrid_pieces.items():
        count = 0
        for kx in range(8):
            for ky in range(8):
                king_pos = (kx, ky)
                
                # Get king's escape squares (adjacent squares)
                escape_squares = set()
                for dx in [-1, 0, 1]:
                    for dy in [-1, 0, 1]:
                        if dx == 0 and dy == 0:
                            continue
                        ex, ey = kx + dx, ky + dy
                        if is_on_board(ex, ey):
                            escape_squares.add((ex, ey))
                
                for hx in range(8):
                    for hy in range(8):
                        hybrid_pos = (hx, hy)
                        if king_pos == hybrid_pos:
                            continue

                        # Get all squares attacked by the hybrid piece
                        attacked_squares = get_piece_attacks(hx, hy)

                        # Condition 1: King must be in check
                        if king_pos not in attacked_squares:
                            continue

                        # Condition 2: All escape squares must be attacked
                        if not escape_squares.issubset(attacked_squares):
                            continue
                            
                        # Condition 3: The hybrid piece cannot be captured by the king
                        if hybrid_pos in escape_squares:
                            continue
                        
                        # If all conditions met, it's a checkmate
                        count += 1
        
        mate_counts[piece_name] = count

    amazon_mates = mate_counts.get("Amazon", 0)
    chancellor_mates = mate_counts.get("Chancellor", 0)
    archbishop_mates = mate_counts.get("Archbishop", 0)
    total_mates = amazon_mates + chancellor_mates + archbishop_mates
    
    print(f"Amazon (Queen + Knight) Mates: {amazon_mates}")
    print(f"Chancellor (Rook + Knight) Mates: {chancellor_mates}")
    print(f"Archbishop (Bishop + Knight) Mates: {archbishop_mates}")
    print(f"Total Distinct Checkmate Positions: {amazon_mates} + {chancellor_mates} + {archbishop_mates} = {total_mates}")

if __name__ == '__main__':
    solve_chess_checkmates()