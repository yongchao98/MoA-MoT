from collections import deque

# Define the board and initial positions
initial_positions = {
    'w1': (1, 1),  # B2
    'w2': (3, 0),  # D1
    'b1': (0, 1),  # A2
    'b2': (2, 0)   # C1
}

# Target positions
target_positions = {
    'w1': (0, 1),  # A2
    'w2': (2, 0),  # C1
    'b1': (1, 1),  # B2
    'b2': (3, 0)   # D1
}

# Possible knight moves (L-shape)
knight_moves = [
    (2, 1), (2, -1), (-2, 1), (-2, -1),
    (1, 2), (1, -2), (-1, 2), (-1, -2)
]

# Convert position to board notation
def pos_to_notation(pos):
    col = chr(pos[0] + ord('A'))
    row = str(pos[1] + 1)
    return col + row

# Check if a position is within the board
def is_within_board(pos):
    return 0 <= pos[0] < 4 and 0 <= pos[1] < 3

# BFS to find the sequence of moves
def find_moves():
    queue = deque([(initial_positions, [])])
    visited = set()
    visited.add(tuple(initial_positions.values()))

    while queue:
        positions, moves = queue.popleft()

        # Check if we reached the target configuration
        if positions == target_positions:
            return moves

        # Determine whose turn it is
        turn = 'w' if len(moves) % 2 == 0 else 'b'

        # Try all possible moves for the current player's knights
        for knight in positions:
            if knight.startswith(turn):
                current_pos = positions[knight]
                for move in knight_moves:
                    new_pos = (current_pos[0] + move[0], current_pos[1] + move[1])
                    if is_within_board(new_pos) and new_pos not in positions.values():
                        new_positions = positions.copy()
                        new_positions[knight] = new_pos
                        new_moves = moves + [f"{turn},{pos_to_notation(current_pos)},{pos_to_notation(new_pos)}"]
                        new_positions_tuple = tuple(new_positions.values())
                        if new_positions_tuple not in visited:
                            visited.add(new_positions_tuple)
                            queue.append((new_positions, new_moves))

    return "No"

# Find the sequence of moves
result = find_moves()
print(result)