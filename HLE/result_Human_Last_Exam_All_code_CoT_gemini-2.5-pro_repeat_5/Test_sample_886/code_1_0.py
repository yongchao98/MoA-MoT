import functools
import collections

# Global constant for board size
BOARD_SIZE = 8

@functools.lru_cache(maxsize=None)
def is_on_board(r, c):
    """Checks if a square (row, col) is on the 8x8 board."""
    return 0 <= r < BOARD_SIZE and 0 <= c < BOARD_SIZE

@functools.lru_cache(maxsize=None)
def get_attacks(piece_char, r, c):
    """
    Gets all squares attacked by a given standard piece from a given square.
    Uses a cache for efficiency.
    """
    attacks = set()
    if piece_char == 'P':  # Pawn (directional attack, assuming a "White" piece moving up the board)
        for dc in [-1, 1]:
            if is_on_board(r + 1, c + dc):
                attacks.add((r + 1, c + dc))
    elif piece_char == 'N':  # Knight
        for dr, dc in [(1, 2), (1, -2), (-1, 2), (-1, -2), (2, 1), (2, -1), (-2, 1), (-2, -1)]:
            if is_on_board(r + dr, c + dc):
                attacks.add((r + dr, c + dc))
    elif piece_char == 'B':  # Bishop
        for dr, dc in [(1, 1), (1, -1), (-1, 1), (-1, -1)]:
            for i in range(1, BOARD_SIZE):
                nr, nc = r + i * dr, c + i * dc
                if is_on_board(nr, nc):
                    attacks.add((nr, nc))
                else:
                    break
    elif piece_char == 'R':  # Rook
        for dr, dc in [(1, 0), (-1, 0), (0, 1), (0, -1)]:
            for i in range(1, BOARD_SIZE):
                nr, nc = r + i * dr, c + i * dc
                if is_on_board(nr, nc):
                    attacks.add((nr, nc))
                else:
                    break
    elif piece_char == 'K':  # King
        for dr in [-1, 0, 1]:
            for dc in [-1, 0, 1]:
                if dr == 0 and dc == 0:
                    continue
                if is_on_board(r + dr, c + dc):
                    attacks.add((r + dr, c + dc))
    return attacks

@functools.lru_cache(maxsize=None)
def get_combined_attacks(piece_combo, r, c):
    """
    Gets all squares attacked by a combined piece (e.g., 'N+B').
    Handles special cases like Queen ('Q') and Amazon ('N+Q').
    """
    total_attacks = set()
    # A Queen ('Q') is a B+R. An Amazon ('N+Q') is N+B+R.
    if 'Q' in piece_combo:
        piece_combo = piece_combo.replace('Q', 'B+R')
        
    components = set(piece_combo.split('+'))
    for component in components:
        total_attacks.update(get_attacks(component, r, c))
    return total_attacks

def solve():
    """
    Calculates the total number of distinct checkmate positions for various custom pieces
    formed by combining the movements of two standard pieces.
    """
    # The list of unique, powerful pieces created by combining two distinct standard pieces.
    # Redundant combinations (e.g., P+K becomes K, B+Q becomes Q) are excluded.
    # 'Q' is included as it represents the B+R combination.
    combined_pieces = [
        ('Pawn+Knight', 'P+N'),
        ('Pawn+Bishop', 'P+B'),
        ('Pawn+Rook', 'P+R'),
        ('Knight+Bishop (Archbishop)', 'N+B'),
        ('Knight+Rook (Chancellor)', 'N+R'),
        ('Knight+Queen (Amazon)', 'N+Q'),
        ('Knight+King', 'N+K'),
        ('Bishop+King', 'B+K'),
        ('Rook+King', 'R+K'),
        ('Queen (Rook+Bishop)', 'Q'),
    ]

    total_mate_positions = 0
    equation_parts = []

    print("Calculating the number of checkmate positions for each custom piece type:")
    
    for name, combo in combined_pieces:
        piece_mate_positions = 0
        # Iterate through all possible king positions (kr, kc)
        for kr in range(BOARD_SIZE):
            for kc in range(BOARD_SIZE):
                king_pos = (kr, kc)
                escape_squares = get_attacks('K', kr, kc)

                # Iterate through all possible piece positions (pr, pc)
                for pr in range(BOARD_SIZE):
                    for pc in range(BOARD_SIZE):
                        piece_pos = (pr, pc)

                        # The piece and king cannot be on the same square.
                        # The piece cannot be captured by the king.
                        if piece_pos == king_pos or piece_pos in escape_squares:
                            continue

                        # Get all squares attacked by the custom piece
                        attacked_by_piece = get_combined_attacks(combo, pr, pc)

                        # Check for mate: King is in check AND all escape squares are covered.
                        if king_pos in attacked_by_piece and escape_squares.issubset(attacked_by_piece):
                            piece_mate_positions += 1
        
        print(f"{name}: {piece_mate_positions}")
        total_mate_positions += piece_mate_positions
        equation_parts.append(str(piece_mate_positions))

    print("\nFinal calculation:")
    print(" + ".join(equation_parts), f"= {total_mate_positions}")
    print(f"\nThe total number of distinct checkmate positions is {total_mate_positions}.")
    
    print(f"<<<{total_mate_positions}>>>")

solve()