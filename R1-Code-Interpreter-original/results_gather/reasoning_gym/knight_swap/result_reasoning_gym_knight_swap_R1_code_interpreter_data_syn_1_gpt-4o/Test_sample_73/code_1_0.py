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
    white_knights = {'C3', 'A2'}
    black_knights = {'B1', 'D1'}
    empty_squares = {'A1', 'B2', 'C1', 'D2', 'D3'}

    # Target positions
    target_white = {'B1', 'D1'}
    target_black = {'C3', 'A2'}

    # Queue for BFS
    queue = deque()
    queue.append((white_knights, black_knights, empty_squares, []))

    # Visited states
    visited = set()
    visited.add((frozenset(white_knights), frozenset(black_knights)))

    while queue:
        w_knights, b_knights, empties, moves = queue.popleft()

        # Check if target is reached
        if w_knights == target_white and b_knights == target_black:
            return moves

        # Try moving white knights
        for w in w_knights:
            for move in get_possible_moves(w, empties):
                new_w_knights = w_knights - {w} | {move}
                new_empties = empties - {move} | {w}
                new_state = (frozenset(new_w_knights), frozenset(b_knights))
                if new_state not in visited:
                    visited.add(new_state)
                    queue.append((new_w_knights, b_knights, new_empties, moves + [f"w,{w},{move}"]))

        # Try moving black knights
        for b in b_knights:
            for move in get_possible_moves(b, empties):
                new_b_knights = b_knights - {b} | {move}
                new_empties = empties - {move} | {b}
                new_state = (frozenset(w_knights), frozenset(new_b_knights))
                if new_state not in visited:
                    visited.add(new_state)
                    queue.append((w_knights, new_b_knights, new_empties, moves + [f"B,{b},{move}"]))

    return "No"

# Execute the function and print the result
print(bfs_knight_swap())