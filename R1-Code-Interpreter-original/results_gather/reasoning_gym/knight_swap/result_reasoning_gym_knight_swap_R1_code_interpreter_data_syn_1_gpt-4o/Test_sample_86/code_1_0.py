from collections import deque

# Define the board and initial positions
board = [
    ['B', 'w', 'w', '.'],
    ['.', '.', '.', '.'],
    ['.', 'B', '.', '.'],
    ['.', '.', '.', '.']
]

# Define the target positions for white and black knights
target_positions = {
    'w': [(2, 0), (1, 1)],  # A3, B1
    'B': [(0, 1), (0, 2)]   # C3, B3
}

# Define knight moves (L-shape)
knight_moves = [
    (2, 1), (2, -1), (-2, 1), (-2, -1),
    (1, 2), (1, -2), (-1, 2), (-1, -2)
]

# Convert board positions to coordinates
def pos_to_coord(pos):
    col, row = pos
    return (3 - (ord(row) - ord('1')), ord(col) - ord('A'))

# Convert coordinates to board positions
def coord_to_pos(coord):
    row, col = coord
    return (chr(col + ord('A')), str(3 - row + 1))

# Check if a position is within the board
def is_within_board(x, y):
    return 0 <= x < 4 and 0 <= y < 4

# BFS to find the sequence of moves
def bfs(start_positions, target_positions):
    queue = deque([(start_positions, [], 'w')])  # (positions, moves, current_turn)
    visited = set()

    while queue:
        positions, moves, turn = queue.popleft()
        if positions == target_positions:
            return moves

        if tuple(positions) in visited:
            continue
        visited.add(tuple(positions))

        for i, (x, y) in enumerate(positions):
            if (turn == 'w' and board[x][y] == 'w') or (turn == 'B' and board[x][y] == 'B'):
                for dx, dy in knight_moves:
                    nx, ny = x + dx, y + dy
                    if is_within_board(nx, ny) and board[nx][ny] == '.':
                        new_positions = positions[:]
                        new_positions[i] = (nx, ny)
                        new_moves = moves + [f"{turn},{coord_to_pos((x, y))},{coord_to_pos((nx, ny))}"]
                        queue.append((new_positions, new_moves, 'B' if turn == 'w' else 'w'))

    return None

# Initial positions of knights
initial_positions = [(0, 0), (0, 1), (2, 1), (1, 1)]  # A3, B3, C3, B1

# Find the sequence of moves
result = bfs(initial_positions, target_positions['w'] + target_positions['B'])

# Output the result
if result:
    print(result)
else:
    print("No")