from collections import deque

# Define the board and initial positions
initial_positions = {
    'w': [(1, 1), (3, 2)],  # B1, D2
    'B': [(1, 3), (3, 1)]   # B3, D1
}

# Target positions
target_positions = {
    'w': [(1, 3), (3, 1)],  # B3, D1
    'B': [(1, 1), (3, 2)]   # B1, D2
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
    col, row = coord
    return (chr(col + ord('A')), str(row + 1))

# Check if a position is within the board
def is_within_board(pos):
    x, y = pos
    return 0 <= x < 4 and 0 <= y < 3

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
            for i, knight_pos in enumerate(positions[color]):
                for move in knight_moves:
                    new_pos = (knight_pos[0] + move[0], knight_pos[1] + move[1])
                    if is_within_board(new_pos) and new_pos not in positions['w'] and new_pos not in positions['B']:
                        new_positions = {k: v[:] for k, v in positions.items()}
                        new_positions[color][i] = new_pos
                        new_moves = moves + [f"{color},{coord_to_pos(knight_pos)},{coord_to_pos(new_pos)}"]
                        state = (tuple(new_positions['w']), tuple(new_positions['B']))
                        if state not in visited:
                            visited.add(state)
                            queue.append((new_positions, new_moves))
    return "No"

# Find the sequence of moves
result = find_moves()
print(result)