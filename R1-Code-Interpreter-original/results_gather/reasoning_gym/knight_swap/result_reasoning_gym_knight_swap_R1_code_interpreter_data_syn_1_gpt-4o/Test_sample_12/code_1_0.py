from collections import deque

def is_valid_move(start, end):
    dx = abs(ord(end[0]) - ord(start[0]))
    dy = abs(int(end[1]) - int(start[1]))
    return (dx == 2 and dy == 1) or (dx == 1 and dy == 2)

def swap_knights_bfs():
    # Initial positions
    initial_white = ('A2', 'B1')
    initial_black = ('C3', 'C1', 'A3')
    initial_empty = ('B2',)

    # Target positions
    target_white = ('C3', 'C1')
    target_black = ('A2', 'B1')

    # Queue for BFS
    queue = deque()
    queue.append((initial_white, initial_black, initial_empty, [], 'w'))

    # Set to track visited states
    visited = set()
    visited.add((initial_white, initial_black, initial_empty))

    while queue:
        white_knights, black_knights, empty_squares, moves, turn = queue.popleft()

        # Check if target configuration is reached
        if set(white_knights) == set(target_white) and set(black_knights) == set(target_black):
            return moves

        # Generate possible moves
        if turn == 'w':
            for i, wk in enumerate(white_knights):
                for es in empty_squares:
                    if is_valid_move(wk, es):
                        new_white = list(white_knights)
                        new_white[i] = es
                        new_empty = list(empty_squares)
                        new_empty.remove(es)
                        new_empty.append(wk)
                        new_state = (tuple(new_white), black_knights, tuple(new_empty))
                        if new_state not in visited:
                            visited.add(new_state)
                            new_moves = moves + [f"w,{wk},{es}"]
                            queue.append((tuple(new_white), black_knights, tuple(new_empty), new_moves, 'B'))
        else:
            for i, bk in enumerate(black_knights):
                for es in empty_squares:
                    if is_valid_move(bk, es):
                        new_black = list(black_knights)
                        new_black[i] = es
                        new_empty = list(empty_squares)
                        new_empty.remove(es)
                        new_empty.append(bk)
                        new_state = (white_knights, tuple(new_black), tuple(new_empty))
                        if new_state not in visited:
                            visited.add(new_state)
                            new_moves = moves + [f"B,{bk},{es}"]
                            queue.append((white_knights, tuple(new_black), tuple(new_empty), new_moves, 'w'))

    return "No"

# Execute the function and print the result
print(swap_knights_bfs())