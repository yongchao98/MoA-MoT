from collections import deque

# Define the initial and target positions
initial_positions = {'w': [(2, 0), (0, 2)], 'B': [(0, 0), (2, 2)]}
target_positions = {'w': [(0, 0), (2, 2)], 'B': [(2, 0), (0, 2)]}

# Define the knight moves
knight_moves = [(2, 1), (2, -1), (-2, 1), (-2, -1), (1, 2), (1, -2), (-1, 2), (-1, -2)]

# Convert board positions to coordinates
def pos_to_coord(pos):
    return (ord(pos[0]) - ord('A'), int(pos[1]) - 1)

# Convert coordinates to board positions
def coord_to_pos(coord):
    return chr(coord[0] + ord('A')) + str(coord[1] + 1)

# Check if a position is within the board
def is_within_board(x, y):
    return 0 <= x < 3 and 0 <= y < 3

# Convert positions dictionary to a hashable type
def positions_to_hashable(positions):
    return tuple((color, tuple(sorted(pos))) for color, pos in positions.items())

# BFS to find the sequence of moves
def bfs_knight_swap():
    queue = deque([(initial_positions, [])])
    visited = set()

    while queue:
        positions, moves = queue.popleft()
        if positions == target_positions:
            return moves

        # Alternate moves between white and black
        current_color = 'w' if len(moves) % 2 == 0 else 'B'
        for i, (x, y) in enumerate(positions[current_color]):
            for dx, dy in knight_moves:
                nx, ny = x + dx, y + dy
                if is_within_board(nx, ny) and (nx, ny) not in positions['w'] and (nx, ny) not in positions['B']:
                    new_positions = {color: list(pos) for color, pos in positions.items()}
                    new_positions[current_color][i] = (nx, ny)
                    new_positions = {color: tuple(pos) for color, pos in new_positions.items()}
                    hashable_positions = positions_to_hashable(new_positions)
                    if hashable_positions not in visited:
                        visited.add(hashable_positions)
                        new_move = f"{current_color},{coord_to_pos((x, y))},{coord_to_pos((nx, ny))}"
                        queue.append((new_positions, moves + [new_move]))

    return "No"

# Find the solution
solution = bfs_knight_swap()
print(solution)