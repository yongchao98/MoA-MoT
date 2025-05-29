from collections import deque

# Define the board and initial positions
initial_positions = {
    'w': [(0, 0), (3, 1)],  # A1, D2
    'B': [(0, 2), (1, 0)]   # A3, B1
}

# Target positions
target_positions = {
    'w': [(0, 2), (1, 0)],  # A3, B1
    'B': [(0, 0), (3, 1)]   # A1, D2
}

# Possible knight moves (L-shape)
knight_moves = [
    (2, 1), (2, -1), (-2, 1), (-2, -1),
    (1, 2), (1, -2), (-1, 2), (-1, -2)
]

# Convert board positions to coordinates
def pos_to_coord(pos):
    col, row = pos
    return (ord(col) - ord('A'), int(row) - 1)

# Convert coordinates to board positions
def coord_to_pos(coord):
    x, y = coord
    return f"{chr(x + ord('A'))}{y + 1}"

# Check if a position is within the board
def is_within_board(x, y):
    return 0 <= x < 4 and 0 <= y < 3

# BFS to find the sequence of moves
def find_moves():
    queue = deque([(initial_positions, [])])
    visited = set()

    while queue:
        positions, moves = queue.popleft()
        if positions == target_positions:
            return moves

        # Determine whose turn it is
        current_color = 'w' if len(moves) % 2 == 0 else 'B'

        for i, (x, y) in enumerate(positions[current_color]):
            for dx, dy in knight_moves:
                nx, ny = x + dx, y + dy
                if is_within_board(nx, ny) and (nx, ny) not in positions['w'] and (nx, ny) not in positions['B']:
                    new_positions = {color: list(pos) for color, pos in positions.items()}
                    new_positions[current_color][i] = (nx, ny)
                    new_positions_tuple = (tuple(new_positions['w']), tuple(new_positions['B']))
                    if new_positions_tuple not in visited:
                        visited.add(new_positions_tuple)
                        new_moves = moves + [f"{current_color},{coord_to_pos((x, y))},{coord_to_pos((nx, ny))}"]
                        queue.append((new_positions, new_moves))

    return "No"

# Find and print the moves
result = find_moves()
print(result)