from collections import deque

# Define the board and initial positions
initial_positions = {
    'w': [(0, 1), (0, 0)],  # A2, A1
    'B': [(3, 2), (3, 0)]   # D3, D1
}

# Target positions
target_positions = {
    'w': [(3, 2), (3, 0)],  # D3, D1
    'B': [(0, 1), (0, 0)]   # A2, A1
}

# Knight moves (L-shape)
knight_moves = [
    (2, 1), (2, -1), (-2, 1), (-2, -1),
    (1, 2), (1, -2), (-1, 2), (-1, -2)
]

# Convert position to board notation
def pos_to_notation(pos):
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
                        new_moves = moves + [f"{color},{pos_to_notation(pos)},{pos_to_notation(new_pos)}"]
                        state = (tuple(new_positions['w']), tuple(new_positions['B']))
                        if state not in visited:
                            visited.add(state)
                            queue.append((new_positions, new_moves))
    return "No"

# Find the solution
solution = bfs_knight_swap()
print(solution)