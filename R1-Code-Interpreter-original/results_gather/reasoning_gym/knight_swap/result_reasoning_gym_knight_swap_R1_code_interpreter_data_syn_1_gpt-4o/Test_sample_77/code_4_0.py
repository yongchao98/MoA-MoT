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
    initial_state = (('B1', 'D2'), ('A3', 'C2'), ('A1', 'C1', 'B3'), ())
    target_state = (('A3', 'C2'), ('B1', 'D2'))

    queue = deque([initial_state])
    visited = set()

    while queue:
        white_knights, black_knights, empty_squares, moves = queue.popleft()

        if (set(white_knights), set(black_knights)) == target_state:
            return list(moves)

        # Generate moves for white knights
        for i, white in enumerate(white_knights):
            for move in get_possible_moves(white, empty_squares):
                new_white_knights = list(white_knights)
                new_white_knights[i] = move
                new_empty_squares = list(empty_squares)
                new_empty_squares.remove(move)
                new_empty_squares.append(white)
                new_moves = moves + (f'w,{white},{move}',)
                new_state = (tuple(new_white_knights), black_knights, tuple(new_empty_squares), new_moves)
                if new_state not in visited:
                    visited.add(new_state)
                    queue.append(new_state)

        # Generate moves for black knights
        for i, black in enumerate(black_knights):
            for move in get_possible_moves(black, empty_squares):
                new_black_knights = list(black_knights)
                new_black_knights[i] = move
                new_empty_squares = list(empty_squares)
                new_empty_squares.remove(move)
                new_empty_squares.append(black)
                new_moves = moves + (f'B,{black},{move}',)
                new_state = (white_knights, tuple(new_black_knights), tuple(new_empty_squares), new_moves)
                if new_state not in visited:
                    visited.add(new_state)
                    queue.append(new_state)

    return "No"

print(bfs_knight_swap())