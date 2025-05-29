from collections import deque

# Define the board and initial positions
initial_positions = {
    'w': [(0, 1), (1, 0)],  # A2, B1
    'B': [(2, 0), (2, 2)],  # C1, C3
}

# Target positions for each knight
target_positions = {
    'w': [(2, 0), (2, 2)],  # C1, C3
    'B': [(0, 1), (1, 0)],  # A2, B1
}

# Possible knight moves (L-shape)
knight_moves = [
    (2, 1), (2, -1), (-2, 1), (-2, -1),
    (1, 2), (1, -2), (-1, 2), (-1, -2)
]

# Convert board positions to coordinates
def pos_to_coord(pos):
    return (ord(pos[0]) - ord('A'), int(pos[1]) - 1)

# Convert coordinates to board positions
def coord_to_pos(coord):
    return chr(coord[0] + ord('A')) + str(coord[1] + 1)

# Check if a position is within the board
def is_within_board(x, y):
    return 0 <= x < 3 and 0 <= y < 3

# BFS to find the sequence of moves
def find_moves():
    queue = deque([(initial_positions, [], 'w')])  # (positions, moves, current_turn)
    visited = set()

    while queue:
        positions, moves, current_turn = queue.popleft()
        if positions == target_positions:
            return moves

        # Generate possible moves for the current turn
        for i, (x, y) in enumerate(positions[current_turn]):
            for dx, dy in knight_moves:
                nx, ny = x + dx, y + dy
                if is_within_board(nx, ny) and (nx, ny) not in positions['w'] and (nx, ny) not in positions['B']:
                    # Create new positions
                    new_positions = {k: v[:] for k, v in positions.items()}
                    new_positions[current_turn][i] = (nx, ny)
                    new_positions_tuple = (tuple(new_positions['w']), tuple(new_positions['B']))

                    if new_positions_tuple not in visited:
                        visited.add(new_positions_tuple)
                        move = f"{current_turn},{coord_to_pos((x, y))},{coord_to_pos((nx, ny))}"
                        queue.append((new_positions, moves + [move], 'B' if current_turn == 'w' else 'w'))

    return "No"

# Find and print the moves
result = find_moves()
print(result)