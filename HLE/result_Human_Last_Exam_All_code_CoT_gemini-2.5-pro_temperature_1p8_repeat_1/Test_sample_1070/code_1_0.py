import collections

def to_coords(alg):
    """Converts algebraic notation (e.g., 'd5') to (file, rank) coordinates."""
    file = ord(alg[0]) - ord('a')
    rank = int(alg[1]) - 1
    return file, rank

def to_alg(pos):
    """Converts (file, rank) coordinates to algebraic notation."""
    file, rank = pos
    return chr(ord('a') + file) + str(rank + 1)

def is_light_square(pos):
    """Checks if a square is light. a1 is dark, so light squares have an odd sum of file and rank."""
    return (pos[0] + pos[1]) % 2 != 0

def get_king_moves(pos):
    """Generates all potential one-square moves for a king."""
    x, y = pos
    for dx in [-1, 0, 1]:
        for dy in [-1, 0, 1]:
            if dx == 0 and dy == 0:
                continue
            nx, ny = x + dx, y + dy
            if 0 <= nx < 8 and 0 <= ny < 8:
                yield (nx, ny)

def get_bishop_moves(pos, occupied_squares):
    """Generates all moves for a bishop from a given position."""
    x, y = pos
    for dx, dy in [(1, 1), (1, -1), (-1, 1), (-1, -1)]:
        for i in range(1, 8):
            nx, ny = x + i * dx, y + i * dy
            if not (0 <= nx < 8 and 0 <= ny < 8):
                break
            yield (nx, ny)
            if (nx, ny) in occupied_squares:
                break

def is_attacked(square, pieces):
    """Checks if a given square is attacked by white."""
    # Attacked by White King
    if max(abs(square[0] - pieces['WK'][0]), abs(square[1] - pieces['WK'][1])) == 1:
        return True
    # Attacked by White Bishop
    wb_pos = pieces['WB']
    occupied = {pieces['WK'], pieces['BK']}
    for move in get_bishop_moves(wb_pos, occupied):
        if move == square:
            return True
    return False

def find_shortest_mate():
    """Performs a BFS to find the shortest mate."""
    start_pieces = {'WK': to_coords('d2'), 'WB': to_coords('e2'), 'BK': to_coords('d5')}
    # A queue of (pieces_tuple, path_list)
    q = collections.deque([(tuple(sorted(start_pieces.items())), [])])
    # A visited set of (pieces_tuple, turn_string)
    visited = {(tuple(sorted(start_pieces.items())), 'W')}

    while q:
        current_pieces_tuple, path = q.popleft()
        current_pieces = dict(current_pieces_tuple)
        turn = 'W' if len(path) % 2 == 0 else 'B'

        if turn == 'W':
            # Try king moves
            for move in get_king_moves(current_pieces['WK']):
                if move == current_pieces['WB']: continue
                if max(abs(move[0] - current_pieces['BK'][0]), abs(move[1] - current_pieces['BK'][1])) <= 1: continue
                
                next_pieces = current_pieces.copy()
                next_pieces['WK'] = move
                
                # Check for mate
                bk_pos = next_pieces['BK']
                if is_attacked(bk_pos, next_pieces):
                    is_mated = True
                    for bk_move in get_king_moves(bk_pos):
                        if is_light_square(bk_move) and not is_attacked(bk_move, next_pieces):
                            is_mated = False
                            break
                    if is_mated:
                        first_move_str = f"K{to_alg(move)}"
                        return first_move_str, (len(path) + 2) // 2

                # Enqueue if not visited
                next_pieces_tuple = tuple(sorted(next_pieces.items()))
                if (next_pieces_tuple, 'B') not in visited:
                    visited.add((next_pieces_tuple, 'B'))
                    q.append((next_pieces_tuple, path + [f"K{to_alg(move)}"]))

            # Try bishop moves
            wb_pos = current_pieces['WB']
            occupied_by_white = {current_pieces['WK']}
            for move in get_bishop_moves(wb_pos, occupied_by_white):
                next_pieces = current_pieces.copy()
                next_pieces['WB'] = move

                # Check for mate
                bk_pos = next_pieces['BK']
                if is_attacked(bk_pos, next_pieces):
                    is_mated = True
                    for bk_move in get_king_moves(bk_pos):
                         if is_light_square(bk_move) and not is_attacked(bk_move, next_pieces):
                            is_mated = False
                            break
                    if is_mated:
                        first_move_str = f"B{to_alg(move)}"
                        return first_move_str, (len(path) + 2) // 2
                
                # Enqueue if not visited
                next_pieces_tuple = tuple(sorted(next_pieces.items()))
                if (next_pieces_tuple, 'B') not in visited:
                    visited.add((next_pieces_tuple, 'B'))
                    q.append((next_pieces_tuple, path + [f"B{to_alg(move)}"]))

        else: # Black's turn
            bk_pos = current_pieces['BK']
            for move in get_king_moves(bk_pos):
                if is_light_square(move) and not is_attacked(move, current_pieces):
                    next_pieces = current_pieces.copy()
                    next_pieces['BK'] = move
                    next_pieces_tuple = tuple(sorted(next_pieces.items()))
                    if (next_pieces_tuple, 'W') not in visited:
                        visited.add((next_pieces_tuple, 'W'))
                        q.append((next_pieces_tuple, path + [f"K{to_alg(move)}"]))
    return None, None

first_move, num_moves = find_shortest_mate()
print(f"{first_move}, {num_moves}")
>>> Kd3, 4