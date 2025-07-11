import collections

def is_on_board(x, y):
    """Checks if a coordinate is within the 8x8 board limits."""
    return 0 <= x < 8 and 0 <= y < 8

def get_pawn_attacks(x, y):
    """Returns the attack squares for a pawn at (x, y).
    Assumes a fixed 'forward' direction (increasing y).
    """
    attacks = set()
    # Diagonal forward-left
    if is_on_board(x - 1, y + 1):
        attacks.add((x - 1, y + 1))
    # Diagonal forward-right
    if is_on_board(x + 1, y + 1):
        attacks.add((x + 1, y + 1))
    return attacks

def get_knight_attacks(x, y):
    """Returns the attack squares for a knight."""
    attacks = set()
    moves = [(-2, -1), (-2, 1), (-1, -2), (-1, 2),
             (1, -2), (1, 2), (2, -1), (2, 1)]
    for dx, dy in moves:
        if is_on_board(x + dx, y + dy):
            attacks.add((x + dx, y + dy))
    return attacks

def get_bishop_attacks(x, y):
    """Returns the attack squares for a bishop."""
    attacks = set()
    for dx, dy in [(-1, -1), (-1, 1), (1, -1), (1, 1)]:
        nx, ny = x + dx, y + dy
        while is_on_board(nx, ny):
            attacks.add((nx, ny))
            nx, ny = nx + dx, ny + dy
    return attacks

def get_rook_attacks(x, y):
    """Returns the attack squares for a rook."""
    attacks = set()
    for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
        nx, ny = x + dx, y + dy
        while is_on_board(nx, ny):
            attacks.add((nx, ny))
            nx, ny = nx + dx, ny + dy
    return attacks

def get_queen_attacks(x, y):
    """Returns the attack squares for a queen."""
    return get_rook_attacks(x, y).union(get_bishop_attacks(x, y))

def get_king_attacks(x, y):
    """Returns the attack/escape squares for a king."""
    attacks = set()
    for dx in [-1, 0, 1]:
        for dy in [-1, 0, 1]:
            if dx == 0 and dy == 0:
                continue
            if is_on_board(x + dx, y + dy):
                attacks.add((x + dx, y + dy))
    return attacks

def solve_chess_mate_puzzle():
    """
    Finds all unassisted checkmate positions for combined-move pieces.
    """
    PIECES = ['P', 'N', 'B', 'R', 'Q', 'K']
    PIECE_ATTACKS = {
        'P': get_pawn_attacks, 'N': get_knight_attacks, 'B': get_bishop_attacks,
        'R': get_rook_attacks, 'Q': get_queen_attacks, 'K': get_king_attacks
    }

    # Generate all 15 combinations of two distinct pieces
    combinations = []
    for i in range(len(PIECES)):
        for j in range(i + 1, len(PIECES)):
            combinations.append((PIECES[i], PIECES[j]))

    total_checkmates = 0
    equation_parts = []

    # Analyze each super-piece combination
    for p1_char, p2_char in combinations:
        combo_mates = 0
        get_p1_attacks = PIECE_ATTACKS[p1_char]
        get_p2_attacks = PIECE_ATTACKS[p2_char]

        # Iterate over all possible king positions (kx, ky)
        for kx in range(8):
            for ky in range(8):
                king_pos = (kx, ky)
                escape_squares = get_king_attacks(kx, ky)
                squares_to_cover = escape_squares.union({king_pos})

                # Iterate over all super-piece positions (spx, spy)
                for spx in range(8):
                    for spy in range(8):
                        sp_pos = (spx, spy)

                        if sp_pos in squares_to_cover:
                            continue

                        # Combine attacks from both movement types
                        sp_attacks = get_p1_attacks(spx, spy).union(get_p2_attacks(spx, spy))

                        # Check for the checkmate condition
                        if king_pos in sp_attacks and escape_squares.issubset(sp_attacks):
                            combo_mates += 1
        
        if combo_mates > 0:
            equation_parts.append(str(combo_mates))
        
        total_checkmates += combo_mates

    # Format the final output as an equation
    equation_string = " + ".join(equation_parts)
    print(f"{equation_string} = {total_checkmates}")

solve_chess_mate_puzzle()