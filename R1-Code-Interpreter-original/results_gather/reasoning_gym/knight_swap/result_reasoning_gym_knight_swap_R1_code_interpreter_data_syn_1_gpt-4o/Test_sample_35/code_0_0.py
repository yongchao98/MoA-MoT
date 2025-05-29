from collections import deque

# Define the board and initial positions
initial_positions = {
    'w': [(0, 2), (2, 0)],  # A3, C1
    'B': [(0, 0), (2, 1)],  # A1, C2
}

# Target positions for the swap
target_positions = {
    'w': [(0, 0), (2, 1)],  # A1, C2
    'B': [(0, 2), (2, 0)],  # A3, C1
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
    return 0 <= x < 3 and 0 <= y < 3

# BFS to find the sequence of moves
def find_moves():
    queue = deque([(initial_positions, [], 'w')])  # (positions, moves, current_turn)
    visited = set()

    while queue:
        positions, moves, turn = queue.popleft()
        positions_tuple = (tuple(positions['w']), tuple(positions['B']))

        if positions_tuple in visited:
            continue
        visited.add(positions_tuple)

        # Check if the current positions match the target positions
        if positions['w'] == target_positions['w'] and positions['B'] == target_positions['B']:
            return moves

        # Generate possible moves for the current turn
        for i, (x, y) in enumerate(positions[turn]):
            for dx, dy in knight_moves:
                nx, ny = x + dx, y + dy
                if is_within_board(nx, ny) and (nx, ny) not in positions['w'] and (nx, ny) not in positions['B']:
                    new_positions = {k: v[:] for k, v in positions.items()}
                    new_positions[turn][i] = (nx, ny)
                    new_moves = moves + [f"{turn},{coord_to_pos((x, y))},{coord_to_pos((nx, ny))}"]
                    next_turn = 'B' if turn == 'w' else 'w'
                    queue.append((new_positions, new_moves, next_turn))

    return None

# Find the sequence of moves
result = find_moves()
if result:
    print(result)
else:
    print("No")