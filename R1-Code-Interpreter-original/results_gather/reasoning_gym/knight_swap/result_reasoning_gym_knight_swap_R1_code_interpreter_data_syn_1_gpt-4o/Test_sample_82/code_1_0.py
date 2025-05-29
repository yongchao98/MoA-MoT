from collections import deque

# Define the board and initial positions
initial_positions = {
    'w': [(1, 2), (0, 1)],  # B3, A2
    'B': [(0, 0), (1, 0)]   # A1, B1
}

# Target positions for the swap
target_positions = {
    'w': [(0, 0), (1, 0)],  # A1, B1
    'B': [(1, 2), (0, 1)]   # B3, A2
}

# Possible knight moves (dx, dy)
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

        # Check if we reached the target configuration
        if positions == target_positions:
            return moves

        # Generate possible moves for the current player
        current_player = 'w' if len(moves) % 2 == 0 else 'B'
        for i, (x, y) in enumerate(positions[current_player]):
            for dx, dy in knight_moves:
                nx, ny = x + dx, y + dy
                if is_within_board(nx, ny) and (nx, ny) not in positions['w'] and (nx, ny) not in positions['B']:
                    # Make the move
                    new_positions = {k: v[:] for k, v in positions.items()}
                    new_positions[current_player][i] = (nx, ny)
                    new_moves = moves + [f"{current_player},{coord_to_pos((x, y))},{coord_to_pos((nx, ny))}"]

                    # Check if this state has been visited
                    state = (tuple(new_positions['w']), tuple(new_positions['B']))
                    if state not in visited:
                        visited.add(state)
                        queue.append((new_positions, new_moves))

    return "No"

# Find the sequence of moves
result = find_moves()
print(result)