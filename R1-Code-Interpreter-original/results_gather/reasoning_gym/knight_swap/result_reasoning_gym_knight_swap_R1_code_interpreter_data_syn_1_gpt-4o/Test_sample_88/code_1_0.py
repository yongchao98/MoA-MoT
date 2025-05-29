from collections import deque

# Define the board and initial positions
initial_positions = {
    'w': [(2, 2), (1, 0)],  # C3, A2
    'B': [(0, 1), (0, 3)]   # B1, D1
}

# Target positions
target_positions = {
    'w': [(0, 1), (0, 3)],  # B1, D1
    'B': [(2, 2), (1, 0)]   # C3, A2
}

# Knight moves (L-shape)
knight_moves = [
    (2, 1), (2, -1), (-2, 1), (-2, -1),
    (1, 2), (1, -2), (-1, 2), (-1, -2)
]

# Convert board positions to coordinates
def pos_to_coord(pos):
    col, row = pos
    return (3 - int(row), ord(col) - ord('A'))

def coord_to_pos(coord):
    row, col = coord
    return f"{chr(col + ord('A'))}{3 - row}"

# Check if a position is within the board
def is_within_board(pos):
    row, col = pos
    return 0 <= row < 3 and 0 <= col < 4

# BFS to find the sequence of moves
def bfs_knight_swap():
    queue = deque([(initial_positions, [])])
    visited = set()

    while queue:
        positions, moves = queue.popleft()
        if positions == target_positions:
            return moves

        # Determine whose turn it is
        turn = 'B' if len(moves) % 2 == 0 else 'w'

        for i, knight_pos in enumerate(positions[turn]):
            for move in knight_moves:
                new_pos = (knight_pos[0] + move[0], knight_pos[1] + move[1])
                if is_within_board(new_pos) and new_pos not in positions['w'] and new_pos not in positions['B']:
                    new_positions = { 'w': positions['w'][:], 'B': positions['B'][:] }
                    new_positions[turn][i] = new_pos
                    new_positions_tuple = (tuple(new_positions['w']), tuple(new_positions['B']))
                    if new_positions_tuple not in visited:
                        visited.add(new_positions_tuple)
                        new_moves = moves + [f"{turn},{coord_to_pos(knight_pos)},{coord_to_pos(new_pos)}"]
                        queue.append((new_positions, new_moves))

    return "No"

# Execute the BFS to find the solution
result = bfs_knight_swap()
print(result)