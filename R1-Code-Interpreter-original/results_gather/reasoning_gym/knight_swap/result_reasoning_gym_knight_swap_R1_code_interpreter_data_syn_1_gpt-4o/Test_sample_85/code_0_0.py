from collections import deque

# Define the board and initial positions
initial_positions = {
    'w': [(0, 2), (2, 2)],  # A3, C3
    'B': [(1, 0), (3, 1)]   # B1, D2
}

# Target positions
target_positions = {
    'w': [(1, 0), (3, 1)],  # B1, D2
    'B': [(0, 2), (2, 2)]   # A3, C3
}

# Knight moves (L-shape)
knight_moves = [
    (2, 1), (2, -1), (-2, 1), (-2, -1),
    (1, 2), (1, -2), (-1, 2), (-1, -2)
]

# Convert position to board coordinates
def pos_to_coord(pos):
    return (ord(pos[0]) - ord('A'), int(pos[1]) - 1)

# Convert board coordinates to position
def coord_to_pos(coord):
    return chr(coord[0] + ord('A')) + str(coord[1] + 1)

# Check if a position is within the board
def is_within_board(x, y):
    return 0 <= x < 4 and 0 <= y < 3

# BFS to find the sequence of moves
def bfs_knight_swap():
    queue = deque([(initial_positions, [])])
    visited = set()

    while queue:
        current_positions, moves = queue.popleft()
        if current_positions == target_positions:
            return moves

        # Alternate moves between white and black
        for color in ['w', 'B']:
            for i, (x, y) in enumerate(current_positions[color]):
                for dx, dy in knight_moves:
                    nx, ny = x + dx, y + dy
                    if is_within_board(nx, ny):
                        # Check if the new position is empty
                        if (nx, ny) not in current_positions['w'] and (nx, ny) not in current_positions['B']:
                            # Create new positions
                            new_positions = {k: v[:] for k, v in current_positions.items()}
                            new_positions[color][i] = (nx, ny)
                            new_positions_tuple = (tuple(new_positions['w']), tuple(new_positions['B']))
                            if new_positions_tuple not in visited:
                                visited.add(new_positions_tuple)
                                new_move = f"{color},{coord_to_pos((x, y))},{coord_to_pos((nx, ny))}"
                                queue.append((new_positions, moves + [new_move]))

    return "No"

# Find the solution
solution = bfs_knight_swap()
print(solution)