from collections import deque

# Define the board and initial positions
initial_positions = {
    'w': [(1, 2), (1, 1)],  # B3, B2
    'B': [(0, 1), (3, 1)]   # A2, D2
}

# Target positions
target_positions = {
    'w': [(0, 1), (3, 1)],  # A2, D2
    'B': [(1, 2), (1, 1)]   # B3, B2
}

# Possible knight moves (dx, dy)
knight_moves = [
    (2, 1), (2, -1), (-2, 1), (-2, -1),
    (1, 2), (1, -2), (-1, 2), (-1, -2)
]

# Convert board position to string
def pos_to_str(pos):
    return chr(pos[0] + ord('A')) + str(pos[1] + 1)

# Check if a position is within the board
def is_within_board(pos):
    return 0 <= pos[0] < 4 and 0 <= pos[1] < 3

# BFS to find the sequence of moves
def bfs_knight_swap():
    queue = deque([(initial_positions, [])])
    visited = set()

    while queue:
        positions, moves = queue.popleft()
        positions_tuple = (tuple(positions['w']), tuple(positions['B']))

        if positions_tuple in visited:
            continue
        visited.add(positions_tuple)

        # Check if target is reached
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
                        new_moves = moves + [f"{color},{pos_to_str(pos)},{pos_to_str(new_pos)}"]
                        queue.append((new_positions, new_moves))

    return "No"

# Find the solution
solution = bfs_knight_swap()
print(solution)