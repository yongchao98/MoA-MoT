from collections import deque

# Define the board and initial positions
initial_positions = {
    'w': [(1, 2), (2, 2)],  # B3, C3
    'B': [(3, 1), (2, 0)]   # D2, C1
}

# Target positions for the swap
target_positions = {
    'w': [(3, 1), (2, 0)],  # D2, C1
    'B': [(1, 2), (2, 2)]   # B3, C3
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
    return 0 <= x < 4 and 0 <= y < 3

# BFS to find the sequence of moves
def find_moves():
    queue = deque([(initial_positions, [], 'w')])  # (positions, moves, current_turn)
    visited = set()

    while queue:
        positions, moves, turn = queue.popleft()
        if positions == target_positions:
            return moves

        # Generate possible moves for the current turn
        for i, (x, y) in enumerate(positions[turn]):
            for dx, dy in knight_moves:
                nx, ny = x + dx, y + dy
                if is_within_board(nx, ny) and (nx, ny) not in positions['w'] and (nx, ny) not in positions['B']:
                    new_positions = { 'w': positions['w'][:], 'B': positions['B'][:] }
                    new_positions[turn][i] = (nx, ny)
                    new_positions_tuple = (tuple(new_positions['w']), tuple(new_positions['B']))
                    if new_positions_tuple not in visited:
                        visited.add(new_positions_tuple)
                        new_moves = moves + [f"{turn},{coord_to_pos((x, y))},{coord_to_pos((nx, ny))}"]
                        queue.append((new_positions, new_moves, 'B' if turn == 'w' else 'w'))

    return "No"

# Find and print the sequence of moves
result = find_moves()
print(result)