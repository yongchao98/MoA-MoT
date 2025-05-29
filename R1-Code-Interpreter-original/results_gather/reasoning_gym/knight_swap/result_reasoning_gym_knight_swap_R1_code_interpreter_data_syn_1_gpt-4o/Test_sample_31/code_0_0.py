from collections import deque

# Define the board and initial positions
initial_positions = {'w': [(0, 2), (1, 1)], 'B': [(0, 1), (3, 0)]}
target_positions = {'w': [(0, 1), (3, 0)], 'B': [(0, 2), (1, 1)]}

# Define knight moves
knight_moves = [(2, 1), (1, 2), (-1, 2), (-2, 1), (-2, -1), (-1, -2), (1, -2), (2, -1)]

# Convert board positions to coordinates
def pos_to_coord(pos):
    col, row = pos
    return ord(col) - ord('A'), int(row) - 1

# Convert coordinates to board positions
def coord_to_pos(coord):
    x, y = coord
    return chr(x + ord('A')) + str(y + 1)

# Check if a position is within the board
def is_within_board(x, y):
    return 0 <= x < 4 and 0 <= y < 3

# BFS to find the sequence of moves
def find_moves():
    queue = deque([(initial_positions, [])])
    visited = set()

    while queue:
        positions, moves = queue.popleft()
        if positions == target_positions:
            return moves

        # Generate next moves
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
                            move = f"{color},{coord_to_pos((x, y))},{coord_to_pos((nx, ny))}"
                            queue.append((new_positions, moves + [move]))

    return "No"

# Find and print the moves
result = find_moves()
print(result)