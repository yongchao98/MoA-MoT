from collections import deque

# Define the initial positions
initial_positions = {
    'w1': 'A1',
    'w2': 'C1',
    'B1': 'B3',
    'B2': 'D2'
}

# Define the target positions
target_positions = {
    'w1': 'B3',
    'w2': 'D2',
    'B1': 'A1',
    'B2': 'C1'
}

# Define the knight moves
knight_moves = [
    (2, 1), (2, -1), (-2, 1), (-2, -1),
    (1, 2), (1, -2), (-1, 2), (-1, -2)
]

# Convert board positions to coordinates
def pos_to_coord(pos):
    col, row = pos
    return ord(col) - ord('A'), int(row) - 1

# Convert coordinates to board positions
def coord_to_pos(coord):
    col, row = coord
    return chr(col + ord('A')) + str(row + 1)

# Check if a position is within the board
def is_within_board(coord):
    col, row = coord
    return 0 <= col < 4 and 0 <= row < 3

# Generate possible moves for a knight
def generate_moves(pos):
    col, row = pos_to_coord(pos)
    moves = []
    for dc, dr in knight_moves:
        new_pos = (col + dc, row + dr)
        if is_within_board(new_pos):
            moves.append(coord_to_pos(new_pos))
    return moves

# BFS to find the sequence of moves
def bfs_solve():
    queue = deque([(initial_positions, [])])
    visited = set()

    while queue:
        positions, moves = queue.popleft()
        if positions == target_positions:
            return moves

        # Alternate moves between white and black
        for knight, pos in positions.items():
            if (len(moves) % 2 == 0 and knight.startswith('w')) or (len(moves) % 2 == 1 and knight.startswith('B')):
                for new_pos in generate_moves(pos):
                    if new_pos not in positions.values():
                        new_positions = positions.copy()
                        new_positions[knight] = new_pos
                        new_moves = moves + [f"{knight[0]},{pos},{new_pos}"]
                        state = tuple(sorted(new_positions.items()))
                        if state not in visited:
                            visited.add(state)
                            queue.append((new_positions, new_moves))
    return "No"

# Find the solution
solution = bfs_solve()
print(solution)