import collections
from itertools import combinations

def is_on_board(x, y):
    """Checks if a square (x, y) is on the 8x8 board."""
    return 0 <= x < 8 and 0 <= y < 8

def get_pawn_w_attacks(x, y):
    """Attacks for a 'White' pawn (moving towards higher rank index)."""
    attacks = set()
    if is_on_board(x - 1, y + 1): attacks.add((x - 1, y + 1))
    if is_on_board(x + 1, y + 1): attacks.add((x + 1, y + 1))
    return attacks

def get_pawn_b_attacks(x, y):
    """Attacks for a 'Black' pawn (moving towards lower rank index)."""
    attacks = set()
    if is_on_board(x - 1, y - 1): attacks.add((x - 1, y - 1))
    if is_on_board(x + 1, y - 1): attacks.add((x + 1, y - 1))
    return attacks

def get_knight_attacks(x, y):
    """Attacks for a Knight."""
    attacks = set()
    moves = [(1, 2), (1, -2), (-1, 2), (-1, -2),
             (2, 1), (2, -1), (-2, 1), (-2, -1)]
    for dx, dy in moves:
        if is_on_board(x + dx, y + dy):
            attacks.add((x + dx, y + dy))
    return attacks

def get_bishop_attacks(x, y):
    """Attacks for a Bishop."""
    attacks = set()
    for dx, dy in [(1, 1), (1, -1), (-1, 1), (-1, -1)]:
        for i in range(1, 8):
            nx, ny = x + i * dx, y + i * dy
            if is_on_board(nx, ny):
                attacks.add((nx, ny))
            else:
                break
    return attacks

def get_rook_attacks(x, y):
    """Attacks for a Rook."""
    attacks = set()
    for dx, dy in [(1, 0), (-1, 0), (0, 1), (0, -1)]:
        for i in range(1, 8):
            nx, ny = x + i * dx, y + i * dy
            if is_on_board(nx, ny):
                attacks.add((nx, ny))
            else:
                break
    return attacks

def get_queen_attacks(x, y):
    """Attacks for a Queen."""
    return get_rook_attacks(x, y).union(get_bishop_attacks(x, y))

def get_king_attacks(x, y):
    """Attacks for a King."""
    attacks = set()
    for dx in [-1, 0, 1]:
        for dy in [-1, 0, 1]:
            if dx == 0 and dy == 0: continue
            if is_on_board(x + dx, y + dy):
                attacks.add((x + dx, y + dy))
    return attacks

def get_super_piece_attacks(piece_combo, x, y):
    """Gets all attacks for a combined piece."""
    all_attacks = set()
    for piece_char in piece_combo:
        all_attacks.update(ATTACK_FUNCTIONS[piece_char](x, y))
    return all_attacks

def is_checkmate(king_pos, piece_pos, piece_combo):
    """Checks if a given board state is checkmate."""
    kx, ky = king_pos
    
    # Get all squares attacked by the super-piece
    piece_attacks = get_super_piece_attacks(piece_combo, piece_pos[0], piece_pos[1])

    # 1. King must be in check
    if king_pos not in piece_attacks:
        return False

    # 2. King must have no legal moves
    for dx in [-1, 0, 1]:
        for dy in [-1, 0, 1]:
            if dx == 0 and dy == 0:
                continue
            
            nkx, nky = kx + dx, ky + dy
            
            # Check if the escape square is on the board
            if not is_on_board(nkx, nky):
                continue
            
            # King cannot move to a square occupied by the attacking piece
            if (nkx, nky) == piece_pos:
                continue

            # King cannot move to a square that is attacked
            if (nkx, nky) not in piece_attacks:
                return False # Found a legal move

    return True # King is in check and has no legal moves

# --- Main Execution ---

# Map piece characters to their attack functions
ATTACK_FUNCTIONS = {
    'P_w': get_pawn_w_attacks,
    'P_b': get_pawn_b_attacks,
    'N': get_knight_attacks,
    'B': get_bishop_attacks,
    'R': get_rook_attacks,
    'Q': get_queen_attacks,
    'K': get_king_attacks,
}

# Define the standard pieces and generate combinations
base_pieces = ['P', 'N', 'B', 'R', 'Q', 'K']
piece_pairs = list(combinations(base_pieces, 2))

# Create the list of super-piece combinations to test
super_piece_combos = []
for p1, p2 in piece_pairs:
    if p1 == 'P':
        super_piece_combos.append(((f'{p1}_w', p2), f"{p1}(w)+{p2}"))
        super_piece_combos.append(((f'{p1}_b', p2), f"{p1}(b)+{p2}"))
    else:
        super_piece_combos.append(((p1, p2), f"{p1}+{p2}"))

# Store results
mate_counts = collections.OrderedDict()
board_squares = [(x, y) for x in range(8) for y in range(8)]

# Iterate through all combinations and board positions
for combo, name in super_piece_combos:
    count = 0
    for king_pos in board_squares:
        for piece_pos in board_squares:
            if king_pos == piece_pos:
                continue
            if is_checkmate(king_pos, piece_pos, combo):
                count += 1
    mate_counts[name] = count

# Print the results in the required format
print("Checkmate positions for each piece combination:")
for name, count in mate_counts.items():
    print(f"- {name}: {count}")

print("\nFinal equation:")
counts_list = list(mate_counts.values())
total_mates = sum(counts_list)
equation_str = " + ".join(map(str, counts_list))
print(f"{equation_str} = {total_mates}")
<<<456>>>