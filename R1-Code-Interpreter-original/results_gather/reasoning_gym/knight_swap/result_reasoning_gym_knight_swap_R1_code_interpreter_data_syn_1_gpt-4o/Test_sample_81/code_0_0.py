from collections import deque

# Define the board and initial positions
initial_positions = {
    'w': [(2, 2), (3, 3)],
    'B': [(2, 3), (3, 2)]
}

# Target positions for the swap
target_positions = {
    'w': [(2, 3), (3, 2)],
    'B': [(2, 2), (3, 3)]
}

# Possible knight moves (L-shape)
knight_moves = [
    (2, 1), (2, -1), (-2, 1), (-2, -1),
    (1, 2), (1, -2), (-1, 2), (-1, -2)
]

# Check if a position is within the board
def is_within_board(x, y):
    return 1 <= x <= 3 and 1 <= y <= 4

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
            for i, (x, y) in enumerate(positions[color]):
                for dx, dy in knight_moves:
                    nx, ny = x + dx, y + dy
                    if is_within_board(nx, ny):
                        # Check if the new position is empty
                        if (nx, ny) not in positions['w'] and (nx, ny) not in positions['B']:
                            # Create new positions
                            new_positions = {k: v[:] for k, v in positions.items()}
                            new_positions[color][i] = (nx, ny)
                            new_positions_tuple = (tuple(new_positions['w']), tuple(new_positions['B']))
                            if new_positions_tuple not in visited:
                                visited.add(new_positions_tuple)
                                new_moves = moves + [f"{color},{chr(64+y)}{x},{chr(64+ny)}{nx}"]
                                queue.append((new_positions, new_moves))
    return "No"

# Find the solution
solution = bfs_knight_swap()
print(solution)