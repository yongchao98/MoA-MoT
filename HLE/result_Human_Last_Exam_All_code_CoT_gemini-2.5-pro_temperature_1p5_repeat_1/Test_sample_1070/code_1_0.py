import collections

def alg_to_coord(alg):
    """Converts algebraic notation like 'd5' to a (row, col) tuple."""
    col = ord(alg[0]) - ord('a')
    row = int(alg[1]) - 1
    return (row, col)

def coord_to_alg(coord):
    """Converts a (row, col) tuple to algebraic notation."""
    row, col = coord
    return f"{chr(ord('a') + col)}{row + 1}"

def is_light_square(coord):
    """Checks if a square is light-colored. a1=(0,0) is dark."""
    row, col = coord
    return (row + col) % 2 != 0

def is_valid(coord):
    """Checks if a coordinate is on the 8x8 board."""
    row, col = coord
    return 0 <= row < 8 and 0 <= col < 8

def get_king_attacks(coord):
    """Returns all squares attacked by a king at a given coordinate."""
    r, c = coord
    moves = set()
    for dr in [-1, 0, 1]:
        for dc in [-1, 0, 1]:
            if dr == 0 and dc == 0:
                continue
            new_coord = (r + dr, c + dc)
            if is_valid(new_coord):
                moves.add(new_coord)
    return moves

def get_bishop_attacks(coord, blockers):
    """Returns all squares attacked by a bishop, considering blockers."""
    r, c = coord
    attacks = set()
    for dr, dc in [(1, 1), (1, -1), (-1, 1), (-1, -1)]:
        for i in range(1, 8):
            new_coord = (r + i * dr, c + i * dc)
            if not is_valid(new_coord):
                break
            attacks.add(new_coord)
            if new_coord in blockers:
                break
    return attacks

def solve_chess_endgame():
    """
    Finds the fastest checkmate for the given chess position and rules.
    """
    # Initial positions
    initial_wk = alg_to_coord('d2')
    initial_wb = alg_to_coord('e2')
    initial_bk = alg_to_coord('d5')

    # State is (wk_pos, wb_pos, bk_pos)
    initial_state = (initial_wk, initial_wb, initial_bk)

    # Queue for BFS: stores (state, path_of_moves)
    queue = collections.deque([(initial_state, [])])
    
    # Visited set to avoid cycles: stores states
    visited = {initial_state}

    while queue:
        (wk_pos, wb_pos, bk_pos), path = queue.popleft()
        
        # It's White's turn if the path length is even
        if len(path) % 2 == 0:
            # Generate White's moves (King)
            for move_coord in get_king_attacks(wk_pos):
                if move_coord not in get_king_attacks(bk_pos) and move_coord != wb_pos:
                    new_state = (move_coord, wb_pos, bk_pos)
                    if new_state not in visited:
                        move_notation = f"K{coord_to_alg(move_coord)}"
                        new_path = path + [move_notation]
                        queue.append((new_state, new_path))
                        visited.add(new_state)

            # Generate White's moves (Bishop)
            blockers = {wk_pos, bk_pos}
            for move_coord in get_bishop_attacks(wb_pos, blockers):
                 # Can't move to an occupied square (except capture, but BK is not captured)
                if move_coord != wk_pos:
                    new_state = (wk_pos, move_coord, bk_pos)
                    if new_state not in visited:
                        move_notation = f"B{coord_to_alg(move_coord)}"
                        new_path = path + [move_notation]
                        queue.append((new_state, new_path))
                        visited.add(new_state)
        else: # Black's turn
            # Check for checkmate first. White's last move might have been mating.
            attacked_by_bishop = get_bishop_attacks(wb_pos, [wk_pos])
            attacked_by_king = get_king_attacks(wk_pos)
            all_attacks = attacked_by_bishop.union(attacked_by_king)

            is_check = bk_pos in attacked_by_bishop

            # Generate Black King's legal moves
            legal_bk_moves = []
            for move_coord in get_king_attacks(bk_pos):
                if is_light_square(move_coord) and move_coord not in all_attacks:
                    legal_bk_moves.append(move_coord)

            if is_check and not legal_bk_moves:
                # Checkmate found
                first_move = path[0]
                num_moves = (len(path) + 1) // 2
                print(f"{first_move}, {num_moves}")
                return

            # If not mate, generate Black's next possible states
            for move_coord in legal_bk_moves:
                new_state = (wk_pos, wb_pos, move_coord)
                if new_state not in visited:
                    move_notation = f"K{coord_to_alg(move_coord)}"
                    new_path = path + [move_notation]
                    queue.append((new_state, new_path))
                    visited.add(new_state)

solve_chess_endgame()