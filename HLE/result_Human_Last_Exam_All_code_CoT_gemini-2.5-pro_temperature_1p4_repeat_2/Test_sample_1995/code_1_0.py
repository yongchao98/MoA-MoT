import collections

def to_coords(s):
    """Converts chess notation like 'a1' to a (row, col) tuple (0-indexed)."""
    col = ord(s[0]) - ord('a')
    row = int(s[1]) - 1
    return row, col

def to_notation(r, c):
    """Converts a (row, col) tuple to chess notation."""
    return f"{chr(ord('a') + c)}{r + 1}"

def get_queen_attacks(r, c, occupied_squares):
    """Gets all squares attacked by a queen at (r, c)."""
    attacks = set()
    # Rook and Bishop directions
    directions = [(1, 0), (-1, 0), (0, 1), (0, -1), (1, 1), (1, -1), (-1, 1), (-1, -1)]
    for dr, dc in directions:
        for i in range(1, 8):
            nr, nc = r + i * dr, c + i * dc
            if 0 <= nr < 8 and 0 <= nc < 8:
                attacks.add((nr, nc))
                if (nr, nc) in occupied_squares:
                    break
            else:
                break
    return attacks

def get_knight_attacks(r, c):
    """Gets all squares attacked by a knight at (r, c)."""
    attacks = set()
    moves = [(1, 2), (1, -2), (-1, 2), (-1, -2),
             (2, 1), (2, -1), (-2, 1), (-2, -1)]
    for dr, dc in moves:
        nr, nc = r + dr, c + dc
        if 0 <= nr < 8 and 0 <= nc < 8:
            attacks.add((nr, nc))
    return attacks

def get_king_neighbors(r, c):
    """Gets all squares adjacent to a king at (r, c)."""
    moves = set()
    for dr in [-1, 0, 1]:
        for dc in [-1, 0, 1]:
            if dr == 0 and dc == 0:
                continue
            nr, nc = r + dr, c + dc
            if 0 <= nr < 8 and 0 <= nc < 8:
                moves.add((nr, nc))
    return moves

def solve_and_verify():
    """
    Verifies the candidate solution for the smallest material to attack 62 squares
    and stalemate a king on the 63rd.
    """
    # Candidate position from chess problem literature
    bk_pos_notation = 'e1'
    wq_pos_notation = 'd3'
    wn_pos_notation = 'e3'

    bk_pos = to_coords(bk_pos_notation)
    wq_pos = to_coords(wq_pos_notation)
    wn_pos = to_coords(wn_pos_notation)

    # For legality, a white king must be on the board. We place it where it
    # does not interfere, e.g., a1. Its material value is not counted.
    # The problem asks for material points of "attacking" pieces.
    wk_pos = to_coords('a1')

    white_pieces = {
        'Q': {'pos': [wq_pos], 'value': 9},
        'N': {'pos': [wn_pos], 'value': 3}
    }

    # All pieces on the board can potentially block attacks.
    occupied = {bk_pos, wq_pos, wn_pos, wk_pos}

    # Calculate all squares attacked by White's pieces
    all_attacked_squares = set()
    for piece_type, data in white_pieces.items():
        for r, c in data['pos']:
            if piece_type == 'Q':
                all_attacked_squares.update(get_queen_attacks(r, c, occupied))
            elif piece_type == 'N':
                all_attacked_squares.update(get_knight_attacks(r, c))

    # --- Verification Step 1: Check for Stalemate ---
    print("Verifying the position...")
    print(f"Black King: {bk_pos_notation}, White Queen: {wq_pos_notation}, White Knight: {wn_pos_notation}\n")

    # 1a. King must not be in check
    if bk_pos in all_attacked_squares:
        print("FAIL: The King is in check. This is not a stalemate position.")
        return

    # 1b. King must have no legal moves
    king_neighbors = get_king_neighbors(bk_pos[0], bk_pos[1])
    unattacked_neighbors = []
    for move in king_neighbors:
        if move not in all_attacked_squares:
            unattacked_neighbors.append(to_notation(move[0], move[1]))

    if unattacked_neighbors:
        print(f"FAIL: The King is not stalemated. It can move to: {', '.join(unattacked_neighbors)}")
        return
    
    print("STALEMATE CONFIRMED:")
    print(f"- The Black King at {bk_pos_notation} is not in check.")
    print(f"- All adjacent squares ({', '.join(sorted([to_notation(r,c) for r,c in king_neighbors]))}) are attacked.\n")


    # --- Verification Step 2: Check Board Coverage ---
    all_squares = set((r, c) for r in range(8) for c in range(8))
    unattacked_by_white = all_squares - all_attacked_squares

    if unattacked_by_white == {bk_pos}:
        print("BOARD COVERAGE CONFIRMED:")
        print(f"- All squares on the board are attacked, except for the King's square ({bk_pos_notation}).\n")
    else:
        other_unattacked = unattacked_by_white - {bk_pos}
        print("FAIL: The board coverage is incorrect.")
        print(f"The following squares are also unattacked: {[to_notation(r,c) for r,c in other_unattacked]}")
        return

    # --- Verification Step 3: Calculate Point Value ---
    print("CALCULATING MATERIAL VALUE:")
    total_value = 0
    calculation_str = []
    for piece_type, data in white_pieces.items():
        value = data['value'] * len(data['pos'])
        piece_name = collections.defaultdict(lambda: piece_type, {'Q':'Queen', 'N':'Knight'})[piece_type]
        print(f"- {piece_name} = {data['value']} points")
        total_value += value
        calculation_str.append(str(data['value']))

    print(f"\nTotal value = {' + '.join(calculation_str)} = {total_value}")


solve_and_verify()