def is_valid_move(start, end):
    dx = abs(ord(end[0]) - ord(start[0]))
    dy = abs(int(end[1]) - int(start[1]))
    return (dx == 2 and dy == 1) or (dx == 1 and dy == 2)

def get_possible_moves(position, empty_squares):
    possible_moves = []
    for dx, dy in [(2, 1), (1, 2), (-1, 2), (-2, 1), (-2, -1), (-1, -2), (1, -2), (2, -1)]:
        new_x = chr(ord(position[0]) + dx)
        new_y = str(int(position[1]) + dy)
        new_position = new_x + new_y
        if new_position in empty_squares:
            possible_moves.append(new_position)
    return possible_moves

def swap_knights(white_knights, black_knights, empty_squares, moves, turn, visited):
    state = (frozenset(white_knights), frozenset(black_knights), frozenset(empty_squares))
    if state in visited:
        return None
    visited.add(state)

    if white_knights == {'A3', 'C1'} and black_knights == {'C3', 'B1'}:
        return moves

    if turn == 'w':
        for knight in list(white_knights):
            possible_moves = get_possible_moves(knight, empty_squares)
            for move in possible_moves:
                white_knights.remove(knight)
                white_knights.add(move)
                empty_squares.remove(move)
                empty_squares.add(knight)
                moves.append(f"w,{knight},{move}")

                result = swap_knights(white_knights, black_knights, empty_squares, moves, 'B', visited)
                if result:
                    return result

                # Backtrack
                moves.pop()
                empty_squares.remove(knight)
                empty_squares.add(move)
                white_knights.remove(move)
                white_knights.add(knight)

    else:
        for knight in list(black_knights):
            possible_moves = get_possible_moves(knight, empty_squares)
            for move in possible_moves:
                black_knights.remove(knight)
                black_knights.add(move)
                empty_squares.remove(move)
                empty_squares.add(knight)
                moves.append(f"B,{knight},{move}")

                result = swap_knights(white_knights, black_knights, empty_squares, moves, 'w', visited)
                if result:
                    return result

                # Backtrack
                moves.pop()
                empty_squares.remove(knight)
                empty_squares.add(move)
                black_knights.remove(move)
                black_knights.add(knight)

    return None

# Initial positions
white_knights = {'C3', 'B1'}
black_knights = {'A3', 'C1'}
empty_squares = {'A1', 'B3', 'D3', 'D2', 'D1'}

# Start the backtracking process
result = swap_knights(white_knights, black_knights, empty_squares, [], 'w', set())
print(result if result else "No")