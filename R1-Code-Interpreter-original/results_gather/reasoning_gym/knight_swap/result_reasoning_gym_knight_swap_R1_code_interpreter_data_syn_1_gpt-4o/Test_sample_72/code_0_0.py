from collections import deque

# Define the board and initial positions
initial_positions = {
    'w': [(3, 3), (2, 0)],  # D3, A2
    'B': [(2, 1), (1, 3)]   # B2, D1
}

# Target positions for the swap
target_positions = {
    'w': [(2, 1), (1, 3)],  # B2, D1
    'B': [(3, 3), (2, 0)]   # D3, A2
}

# Possible knight moves (L-shape)
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
    return chr(col + ord('A')) + str(3 - row)

# Check if a position is within the board
def is_within_board(pos):
    row, col = pos
    return 0 <= row < 4 and 0 <= col < 4

# BFS to find the sequence of moves
def find_moves():
    queue = deque([(initial_positions, [])])
    visited = set()

    while queue:
        positions, moves = queue.popleft()
        if positions == target_positions:
            return moves

        # Alternate moves between white and black
        for color in ['w', 'B']:
            for i, pos in enumerate(positions[color]):
                for move in knight_moves:
                    new_pos = (pos[0] + move[0], pos[1] + move[1])
                    if is_within_board(new_pos) and new_pos not in positions['w'] and new_pos not in positions['B']:
                        new_positions = {k: v[:] for k, v in positions.items()}
                        new_positions[color][i] = new_pos
                        new_moves = moves + [f"{color},{coord_to_pos(pos)},{coord_to_pos(new_pos)}"]
                        state = (tuple(new_positions['w']), tuple(new_positions['B']))
                        if state not in visited:
                            visited.add(state)
                            queue.append((new_positions, new_moves))
    return "No"

# Find and print the moves
result = find_moves()
print(result)