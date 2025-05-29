from collections import deque

# Define the board and initial positions
initial_positions = {
    'w': [(0, 1), (2, 0)],  # A2, C1
    'B': [(3, 2), (3, 1)]   # D3, D2
}

# Target positions for the swap
target_positions = {
    'w': [(3, 2), (3, 1)],  # D3, D2
    'B': [(0, 1), (2, 0)]   # A2, C1
}

# Possible knight moves (L-shape)
knight_moves = [
    (2, 1), (2, -1), (-2, 1), (-2, -1),
    (1, 2), (1, -2), (-1, 2), (-1, -2)
]

# Convert board positions to coordinates
def pos_to_coord(pos):
    return ord(pos[0]) - ord('A'), int(pos[1]) - 1

# Convert coordinates to board positions
def coord_to_pos(coord):
    return chr(coord[0] + ord('A')) + str(coord[1] + 1)

# Check if a position is within the board
def is_within_board(x, y):
    return 0 <= x < 4 and 0 <= y < 3

# BFS to find the sequence of moves
def find_moves():
    queue = deque([(initial_positions, [])])  # (current positions, move sequence)
    visited = set()

    while queue:
        current_positions, moves = queue.popleft()

        # Check if the current positions match the target positions
        if current_positions == target_positions:
            return moves

        # Alternate moves between white and black knights
        for color in ['w', 'B']:
            for i, (x, y) in enumerate(current_positions[color]):
                for dx, dy in knight_moves:
                    nx, ny = x + dx, y + dy
                    if is_within_board(nx, ny):
                        # Check if the new position is empty
                        if (nx, ny) not in current_positions['w'] and (nx, ny) not in current_positions['B']:
                            # Create new positions
                            new_positions = {k: v[:] for k, v in current_positions.items()}
                            new_positions[color][i] = (nx, ny)
                            new_positions_tuple = (tuple(new_positions['w']), tuple(new_positions['B']))

                            if new_positions_tuple not in visited:
                                visited.add(new_positions_tuple)
                                new_move = f"{color},{coord_to_pos((x, y))},{coord_to_pos((nx, ny))}"
                                queue.append((new_positions, moves + [new_move]))

    return "No"

# Find and print the sequence of moves
result = find_moves()
print(result)