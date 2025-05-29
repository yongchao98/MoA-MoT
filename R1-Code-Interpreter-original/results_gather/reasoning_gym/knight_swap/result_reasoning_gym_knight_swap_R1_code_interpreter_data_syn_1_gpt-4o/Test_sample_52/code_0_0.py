from collections import deque

# Define the board and initial positions
initial_positions = {
    'w1': (2, 2),  # C3
    'w2': (1, 1),  # B1
    'b1': (0, 2),  # A3
    'b2': (3, 0)   # D1
}

# Target positions
target_positions = {
    'w1': (0, 2),  # A3
    'w2': (3, 0),  # D1
    'b1': (2, 2),  # C3
    'b2': (1, 1)   # B1
}

# Possible knight moves
knight_moves = [
    (2, 1), (2, -1), (-2, 1), (-2, -1),
    (1, 2), (1, -2), (-1, 2), (-1, -2)
]

# Check if a position is within the board
def is_within_board(x, y):
    return 0 <= x < 4 and 0 <= y < 3

# BFS to find the sequence of moves
def bfs_knight_swap():
    queue = deque([(initial_positions, [])])
    visited = set()
    visited.add(tuple(initial_positions.values()))

    while queue:
        current_positions, moves = queue.popleft()

        # Check if current positions match target positions
        if current_positions == target_positions:
            return moves

        # Determine whose turn it is
        turn = 'w' if len(moves) % 2 == 0 else 'b'

        # Try all possible moves for the current player's knights
        for knight, pos in current_positions.items():
            if knight.startswith(turn):
                x, y = pos
                for dx, dy in knight_moves:
                    new_x, new_y = x + dx, y + dy
                    if is_within_board(new_x, new_y) and (new_x, new_y) not in current_positions.values():
                        new_positions = current_positions.copy()
                        new_positions[knight] = (new_x, new_y)
                        new_positions_tuple = tuple(new_positions.values())
                        if new_positions_tuple not in visited:
                            visited.add(new_positions_tuple)
                            move_str = f"{turn},{chr(65 + x)}{3 - y},{chr(65 + new_x)}{3 - new_y}"
                            queue.append((new_positions, moves + [move_str]))

    return "No"

# Execute the BFS to find the solution
result = bfs_knight_swap()
print(result)