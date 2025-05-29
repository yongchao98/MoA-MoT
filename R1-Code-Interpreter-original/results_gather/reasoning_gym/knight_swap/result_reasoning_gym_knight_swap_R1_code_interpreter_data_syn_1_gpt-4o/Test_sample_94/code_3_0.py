from collections import deque

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

def swap_knights_bfs():
    # Initial positions
    white_knights = {'C3', 'B1'}
    black_knights = {'A3', 'C1'}
    empty_squares = {'A1', 'B3', 'D3', 'D2', 'D1'}

    # Target positions
    target_white = {'A3', 'C1'}
    target_black = {'C3', 'B1'}

    # Queue for BFS
    queue = deque()
    queue.append((white_knights, black_knights, empty_squares, [], 'w'))

    # Visited states
    visited = set()

    while queue:
        white_knights, black_knights, empty_squares, moves, turn = queue.popleft()
        state = (frozenset(white_knights), frozenset(black_knights), frozenset(empty_squares))
        
        if state in visited:
            continue
        visited.add(state)

        if white_knights == target_white and black_knights == target_black:
            return moves

        if turn == 'w':
            for knight in list(white_knights):
                possible_moves = get_possible_moves(knight, empty_squares)
                for move in possible_moves:
                    new_white_knights = white_knights.copy()
                    new_white_knights.remove(knight)
                    new_white_knights.add(move)
                    new_empty_squares = empty_squares.copy()
                    new_empty_squares.remove(move)
                    new_empty_squares.add(knight)
                    new_moves = moves + [f"w,{knight},{move}"]
                    queue.append((new_white_knights, black_knights, new_empty_squares, new_moves, 'B'))
        else:
            for knight in list(black_knights):
                possible_moves = get_possible_moves(knight, empty_squares)
                for move in possible_moves:
                    new_black_knights = black_knights.copy()
                    new_black_knights.remove(knight)
                    new_black_knights.add(move)
                    new_empty_squares = empty_squares.copy()
                    new_empty_squares.remove(move)
                    new_empty_squares.add(knight)
                    new_moves = moves + [f"B,{knight},{move}"]
                    queue.append((white_knights, new_black_knights, new_empty_squares, new_moves, 'w'))

    return "No"

result = swap_knights_bfs()
print(result)