from collections import deque

# Define the board and initial positions
initial_positions = {
    'w': [(1, 2), (2, 2)],  # B3, C3
    'B': [(0, 2), (1, 0)]   # A3, B1
}

# Target positions for the swap
target_positions = {
    'w': [(0, 2), (1, 0)],  # A3, B1
    'B': [(1, 2), (2, 2)]   # B3, C3
}

# Possible knight moves (L-shape)
knight_moves = [
    (2, 1), (2, -1), (-2, 1), (-2, -1),
    (1, 2), (1, -2), (-1, 2), (-1, -2)
]

# Check if a position is within the board
def is_within_board(x, y):
    return 0 <= x < 4 and 0 <= y < 3

# BFS to find the sequence of moves
def find_moves():
    queue = deque([(initial_positions, [])])  # (positions, moves)
    visited = set()

    while queue:
        positions, moves = queue.popleft()
        current_positions = (tuple(positions['w']), tuple(positions['B']))

        if current_positions in visited:
            continue
        visited.add(current_positions)

        # Check if the current positions match the target positions
        if positions['w'] == target_positions['w'] and positions['B'] == target_positions['B']:
            return moves

        # Alternate moves between white and black knights
        for color in ['w', 'B']:
            for i, (x, y) in enumerate(positions[color]):
                for dx, dy in knight_moves:
                    nx, ny = x + dx, y + dy
                    if is_within_board(nx, ny) and (nx, ny) not in positions['w'] and (nx, ny) not in positions['B']:
                        new_positions = {k: v[:] for k, v in positions.items()}
                        new_positions[color][i] = (nx, ny)
                        new_moves = moves + [f"{color},{chr(65+x)}{3-y},{chr(65+nx)}{3-ny}"]
                        queue.append((new_positions, new_moves))

    return "No"

# Find and print the sequence of moves
result = find_moves()
print(result)