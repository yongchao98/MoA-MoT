from collections import deque

# Define the board and initial positions
initial_positions = {
    'w': [(2, 2), (3, 1)],  # C2, D1
    'B': [(1, 1), (4, 2)]   # A1, D2
}

# Target positions
target_positions = {
    'w': [(1, 1), (4, 2)],  # A1, D2
    'B': [(2, 2), (3, 1)]   # C2, D1
}

# Possible knight moves (L-shape)
knight_moves = [
    (2, 1), (2, -1), (-2, 1), (-2, -1),
    (1, 2), (1, -2), (-1, 2), (-1, -2)
]

# Convert board positions to coordinates
def pos_to_coord(pos):
    col, row = pos
    return (ord(col) - ord('A') + 1, int(row))

# Convert coordinates to board positions
def coord_to_pos(coord):
    x, y = coord
    return f"{chr(x + ord('A') - 1)}{y}"

# Check if a position is within the board
def is_within_board(x, y):
    return 1 <= x <= 4 and 1 <= y <= 3

# BFS to find the sequence of moves
def find_moves():
    queue = deque([(initial_positions, [])])
    visited = set()

    while queue:
        positions, moves = queue.popleft()
        current_positions = (tuple(positions['w']), tuple(positions['B']))

        if current_positions in visited:
            continue
        visited.add(current_positions)

        # Check if the current positions match the target positions
        if positions == target_positions:
            return moves

        # Determine whose turn it is (w or B)
        turn = 'w' if len(moves) % 2 == 0 else 'B'

        # Try all possible moves for the current player's knights
        for i, (x, y) in enumerate(positions[turn]):
            for dx, dy in knight_moves:
                nx, ny = x + dx, y + dy
                if is_within_board(nx, ny) and (nx, ny) not in positions['w'] and (nx, ny) not in positions['B']:
                    # Make the move
                    new_positions = { 'w': positions['w'][:], 'B': positions['B'][:] }
                    new_positions[turn][i] = (nx, ny)
                    move = f"{turn},{coord_to_pos((x, y))},{coord_to_pos((nx, ny))}"
                    queue.append((new_positions, moves + [move]))

    return "No"

# Find and print the moves
result = find_moves()
print(result)