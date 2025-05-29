from collections import deque

# Define the board and initial positions
initial_positions = {
    'w': [(1, 2), (0, 1)],  # B3, A2
    'B': [(3, 1), (0, 0)],  # D2, A1
}

# Target positions for the swap
target_positions = {
    'w': [(3, 1), (0, 0)],  # D2, A1
    'B': [(1, 2), (0, 1)],  # B3, A2
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
    x, y = coord
    return (chr(x + ord('A')), y + 1)

# Check if a position is within the board
def is_within_board(x, y):
    return 0 <= x < 4 and 0 <= y < 3

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
                    if is_within_board(nx, ny) and (nx, ny) not in positions['w'] and (nx, ny) not in positions['B']:
                        new_positions = {k: v[:] for k, v in positions.items()}
                        new_positions[color][i] = (nx, ny)
                        new_positions_tuple = (tuple(new_positions['w']), tuple(new_positions['B']))
                        if new_positions_tuple not in visited:
                            visited.add(new_positions_tuple)
                            new_moves = moves + [f"{color},{coord_to_pos((x, y))},{coord_to_pos((nx, ny))}"]
                            queue.append((new_positions, new_moves))
    return "No"

# Find the solution
solution = bfs_knight_swap()
print(solution)