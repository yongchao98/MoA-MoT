from collections import deque

# Define the board and initial positions
board = [
    ['.', 'B', 'w', '.'],
    [' ', 'w', ' ', 'B'],
    [' ', '.', '.', ' ']
]

# Define the target positions for white and black knights
target_white = [(1, 1), (0, 3)]
target_black = [(0, 2), (1, 0)]

# Define knight moves (L-shape)
knight_moves = [
    (2, 1), (2, -1), (-2, 1), (-2, -1),
    (1, 2), (1, -2), (-1, 2), (-1, -2)
]

# Convert board positions to coordinates
def pos_to_coord(pos):
    col, row = pos
    return (3 - int(row), ord(col) - ord('A'))

# Convert coordinates to board positions
def coord_to_pos(coord):
    row, col = coord
    return f"{chr(col + ord('A'))}{3 - row}"

# Find initial positions of knights
white_knights = [pos_to_coord('C3'), pos_to_coord('B2')]
black_knights = [pos_to_coord('B3'), pos_to_coord('D2')]

# BFS to find the sequence of moves
def bfs(white_knights, black_knights):
    queue = deque([(tuple(white_knights), tuple(black_knights), [], 'w')])
    visited = set()
    visited.add((tuple(white_knights), tuple(black_knights)))

    while queue:
        w_pos, b_pos, moves, turn = queue.popleft()

        # Check if the current positions match the target
        if sorted(w_pos) == sorted(target_white) and sorted(b_pos) == sorted(target_black):
            return moves

        # Determine the current player's knights
        current_knights = w_pos if turn == 'w' else b_pos

        # Try all possible moves for the current player's knights
        for i, (r, c) in enumerate(current_knights):
            for dr, dc in knight_moves:
                nr, nc = r + dr, c + dc
                if 0 <= nr < 3 and 0 <= nc < 4 and (nr, nc) not in w_pos and (nr, nc) not in b_pos:
                    new_knights = list(current_knights)
                    new_knights[i] = (nr, nc)
                    new_knights = tuple(sorted(new_knights))
                    new_w_pos = new_knights if turn == 'w' else w_pos
                    new_b_pos = new_knights if turn == 'b' else b_pos
                    if (new_w_pos, new_b_pos) not in visited:
                        visited.add((new_w_pos, new_b_pos))
                        move = f"{turn},{coord_to_pos((r, c))},{coord_to_pos((nr, nc))}"
                        queue.append((new_w_pos, new_b_pos, moves + [move], 'b' if turn == 'w' else 'w'))

    return None

# Find the sequence of moves
result = bfs(white_knights, black_knights)

# Output the result
if result:
    print(result)
else:
    print("No")