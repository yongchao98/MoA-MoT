from collections import deque

# Define the initial and target positions
initial_positions = {'w': [(3, 4), (1, 2)], 'B': [(3, 2), (1, 1)]}
target_positions = {'w': [(3, 2), (1, 1)], 'B': [(3, 4), (1, 2)]}

# Define the knight moves (L-shape)
knight_moves = [(2, 1), (2, -1), (-2, 1), (-2, -1), (1, 2), (1, -2), (-1, 2), (-1, -2)]

# Convert board positions to coordinates
def pos_to_coord(pos):
    col, row = pos
    return (int(row), ord(col) - ord('A') + 1)

# Convert coordinates to board positions
def coord_to_pos(coord):
    row, col = coord
    return (chr(col + ord('A') - 1), str(row))

# Check if a position is within the board
def is_within_board(pos):
    row, col = pos
    return 1 <= row <= 3 and 1 <= col <= 4

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

# Execute the BFS to find the solution
result = bfs_knight_swap()
print(result)