import collections

def get_rook_attacks(pos):
    """Returns all squares attacked by a rook from pos."""
    x, y = pos
    attacks = set()
    for i in range(8):
        if i != x: attacks.add((i, y))
        if i != y: attacks.add((x, i))
    return attacks

def get_bishop_attacks(pos):
    """Returns all squares attacked by a bishop from pos."""
    x, y = pos
    attacks = set()
    for i in range(1, 8):
        # Top-right
        if 0 <= x + i < 8 and 0 <= y + i < 8: attacks.add((x + i, y + i))
        # Top-left
        if 0 <= x - i < 8 and 0 <= y + i < 8: attacks.add((x - i, y + i))
        # Bottom-right
        if 0 <= x + i < 8 and 0 <= y - i < 8: attacks.add((x + i, y - i))
        # Bottom-left
        if 0 <= x - i < 8 and 0 <= y - i < 8: attacks.add((x - i, y - i))
    return attacks

def get_queen_attacks(pos):
    """Returns all squares attacked by a queen from pos."""
    return get_rook_attacks(pos).union(get_bishop_attacks(pos))

def get_knight_attacks(pos):
    """Returns all squares attacked by a knight from pos."""
    x, y = pos
    attacks = set()
    moves = [(1, 2), (1, -2), (-1, 2), (-1, -2),
             (2, 1), (2, -1), (-2, 1), (-2, -1)]
    for dx, dy in moves:
        if 0 <= x + dx < 8 and 0 <= y + dy < 8:
            attacks.add((x + dx, y + dy))
    return attacks

def get_king_attacks(pos):
    """Returns all squares a king can move to from pos."""
    x, y = pos
    attacks = set()
    for dx in [-1, 0, 1]:
        for dy in [-1, 0, 1]:
            if dx == 0 and dy == 0: continue
            if 0 <= x + dx < 8 and 0 <= y + dy < 8:
                attacks.add((x + dx, y + dy))
    return attacks

def get_pawn_attacks(pos):
    """Returns squares attacked by a pawn, assuming forward is increasing y."""
    x, y = pos
    attacks = set()
    # Diagonal forward captures
    if 0 <= x - 1 < 8 and 0 <= y + 1 < 8: attacks.add((x - 1, y + 1))
    if 0 <= x + 1 < 8 and 0 <= y + 1 < 8: attacks.add((x + 1, y + 1))
    return attacks

ATTACK_FUNCTIONS = {
    'K': get_king_attacks,
    'Q': get_queen_attacks,
    'R': get_rook_attacks,
    'B': get_bishop_attacks,
    'N': get_knight_attacks,
    'P': get_pawn_attacks
}

def is_checkmate(combo, h_pos, k_pos):
    """Checks if a hybrid piece at h_pos checkmates a king at k_pos."""
    # Piece cannot be adjacent to the king, otherwise king can capture.
    if max(abs(h_pos[0] - k_pos[0]), abs(h_pos[1] - k_pos[1])) <= 1:
        return False

    # Get the combined attacks of the hybrid piece.
    p1, p2 = combo
    total_attacks = ATTACK_FUNCTIONS[p1](h_pos).union(ATTACK_FUNCTIONS[p2](h_pos))

    # Get the squares the king needs to be prevented from occupying.
    squares_to_control = get_king_attacks(k_pos)
    squares_to_control.add(k_pos)

    # Check if the hybrid piece attacks the king and all escape squares.
    return squares_to_control.issubset(total_attacks)

def solve():
    """
    Calculates the number of checkmate positions for all 2-piece combinations.
    """
    pieces = ['K', 'Q', 'R', 'B', 'N', 'P']
    combinations = []
    for i in range(len(pieces)):
        for j in range(i + 1, len(pieces)):
            combinations.append((pieces[i], pieces[j]))

    board_squares = [(x, y) for x in range(8) for y in range(8)]
    results = collections.OrderedDict()
    total_checkmates = 0

    for combo in combinations:
        combo_name = f"{combo[0]}+{combo[1]}"
        count = 0
        for h_pos in board_squares:
            for k_pos in board_squares:
                if h_pos == k_pos:
                    continue
                if is_checkmate(combo, h_pos, k_pos):
                    count += 1
        results[combo_name] = count
        total_checkmates += count
    
    # Format the final equation string
    equation_parts = []
    for name, count in results.items():
        equation_parts.append(f"{count}")
    
    print("The number of distinct checkmate positions for each combination are:")
    # Print the itemized list
    for name, count in results.items():
        print(f"- {name}: {count}")
    
    # Print the final equation
    final_equation = " + ".join(equation_parts)
    print("\nThe final calculation is:")
    print(f"{final_equation} = {total_checkmates}")

if __name__ == '__main__':
    solve()
    print("\n<<<984>>>")
