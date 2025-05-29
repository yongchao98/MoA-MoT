from collections import deque

def is_valid_move(start, end):
    dx = abs(ord(end[0]) - ord(start[0]))
    dy = abs(int(end[1]) - int(start[1]))
    return (dx == 2 and dy == 1) or (dx == 1 and dy == 2)

def get_possible_moves(position, empty_squares):
    possible_moves = []
    for empty in empty_squares:
        if is_valid_move(position, empty):
            possible_moves.append(empty)
    return possible_moves

def bfs_knight_swap():
    # Initial positions
    initial_white = {'A1': 'w', 'C1': 'w'}
    initial_black = {'B2': 'B', 'D2': 'B'}
    initial_empty = {'A3', 'B3', 'C3', 'D3'}

    # Target positions
    target_white = {'B2': 'w', 'D2': 'w'}
    target_black = {'A1': 'B', 'C1': 'B'}

    # Queue for BFS: (white_positions, black_positions, empty_squares, moves, turn)
    queue = deque([(initial_white, initial_black, initial_empty, [], 'B')])
    visited = set()

    while queue:
        white_knights, black_knights, empty_squares, moves, turn = queue.popleft()

        # Check if target is reached
        if white_knights == target_white and black_knights == target_black:
            return moves

        # Generate next states
        if turn == 'B':
            for pos in black_knights:
                for move in get_possible_moves(pos, empty_squares):
                    new_black_knights = black_knights.copy()
                    new_black_knights[move] = new_black_knights.pop(pos)
                    new_empty_squares = empty_squares.copy()
                    new_empty_squares.remove(move)
                    new_empty_squares.add(pos)
                    new_moves = moves + [f'B,{pos},{move}']
                    state = (frozenset(white_knights.items()), frozenset(new_black_knights.items()), frozenset(new_empty_squares))
                    if state not in visited:
                        visited.add(state)
                        queue.append((white_knights, new_black_knights, new_empty_squares, new_moves, 'w'))
        else:
            for pos in white_knights:
                for move in get_possible_moves(pos, empty_squares):
                    new_white_knights = white_knights.copy()
                    new_white_knights[move] = new_white_knights.pop(pos)
                    new_empty_squares = empty_squares.copy()
                    new_empty_squares.remove(move)
                    new_empty_squares.add(pos)
                    new_moves = moves + [f'w,{pos},{move}']
                    state = (frozenset(new_white_knights.items()), frozenset(black_knights.items()), frozenset(new_empty_squares))
                    if state not in visited:
                        visited.add(state)
                        queue.append((new_white_knights, black_knights, new_empty_squares, new_moves, 'B'))

    return "No"

# Execute the function and print the result
print(bfs_knight_swap())