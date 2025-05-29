from collections import deque

def is_valid_move(start, end):
    dx = abs(ord(end[0]) - ord(start[0]))
    dy = abs(int(end[1]) - int(start[1]))
    return (dx == 2 and dy == 1) or (dx == 1 and dy == 2)

def get_possible_moves(knight, empty_squares):
    possible_moves = []
    for empty in empty_squares:
        if is_valid_move(knight, empty):
            possible_moves.append(empty)
    return possible_moves

def bfs_solve():
    # Initial positions
    initial_white = ['A3', 'D1']
    initial_black = ['B1', 'C1']
    initial_empty = ['B3', 'A1']

    # Target positions
    target_white = ['B1', 'C1']
    target_black = ['A3', 'D1']

    # Queue for BFS
    queue = deque()
    queue.append((initial_white, initial_black, initial_empty, []))

    # Set to track visited states
    visited = set()
    visited.add((tuple(initial_white), tuple(initial_black)))

    while queue:
        white_knights, black_knights, empty_squares, moves = queue.popleft()

        # Check if target configuration is reached
        if set(white_knights) == set(target_white) and set(black_knights) == set(target_black):
            return moves

        # Generate possible moves for white knights
        for i, white in enumerate(white_knights):
            for move in get_possible_moves(white, empty_squares):
                new_white_knights = white_knights[:]
                new_white_knights[i] = move
                new_empty_squares = empty_squares[:]
                new_empty_squares.remove(move)
                new_empty_squares.append(white)
                new_moves = moves + [f"w,{white},{move}"]

                state = (tuple(new_white_knights), tuple(black_knights))
                if state not in visited:
                    visited.add(state)
                    queue.append((new_white_knights, black_knights, new_empty_squares, new_moves))

        # Generate possible moves for black knights
        for i, black in enumerate(black_knights):
            for move in get_possible_moves(black, empty_squares):
                new_black_knights = black_knights[:]
                new_black_knights[i] = move
                new_empty_squares = empty_squares[:]
                new_empty_squares.remove(move)
                new_empty_squares.append(black)
                new_moves = moves + [f"B,{black},{move}"]

                state = (tuple(white_knights), tuple(new_black_knights))
                if state not in visited:
                    visited.add(state)
                    queue.append((white_knights, new_black_knights, new_empty_squares, new_moves))

    return "No"

# Execute the function and print the result
result = bfs_solve()
print(result)