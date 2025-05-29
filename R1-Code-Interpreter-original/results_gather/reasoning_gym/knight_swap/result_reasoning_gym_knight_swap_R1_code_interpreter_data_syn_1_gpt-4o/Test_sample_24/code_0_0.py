from collections import deque

# Define the initial positions of the knights
initial_positions = {
    'w': [(0, 2), (3, 1)],  # A3, D2
    'B': [(1, 2), (2, 1)]   # B3, C2
}

# Define the target positions for the swap
target_positions = {
    'w': [(1, 2), (2, 1)],  # B3, C2
    'B': [(0, 2), (3, 1)]   # A3, D2
}

# Define the possible knight moves
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

# Perform a BFS to find the sequence of moves
def find_moves():
    queue = deque([(initial_positions, [])])
    visited = set()

    while queue:
        positions, moves = queue.popleft()
        if positions == target_positions:
            return moves

        # Generate possible moves for each knight
        for color in ['w', 'B']:
            for i, (x, y) in enumerate(positions[color]):
                for dx, dy in knight_moves:
                    nx, ny = x + dx, y + dy
                    if is_within_board(nx, ny):
                        # Check if the target square is empty
                        if (nx, ny) not in positions['w'] and (nx, ny) not in positions['B']:
                            # Create new positions
                            new_positions = {k: v[:] for k, v in positions.items()}
                            new_positions[color][i] = (nx, ny)
                            new_positions_tuple = (tuple(new_positions['w']), tuple(new_positions['B']))
                            if new_positions_tuple not in visited:
                                visited.add(new_positions_tuple)
                                new_move = f"{color},{coord_to_pos((x, y))},{coord_to_pos((nx, ny))}"
                                queue.append((new_positions, moves + [new_move]))

    return "No"

# Find the sequence of moves
result = find_moves()
print(result)