from collections import deque

# Define the board and initial positions
initial_positions = {
    'w': [(2, 2), (0, 1)],  # C3, A2
    'B': [(0, 0), (2, 0)]   # A1, C1
}

# Target positions for each knight
target_positions = {
    'w': [(0, 0), (2, 0)],  # A1, C1
    'B': [(2, 2), (0, 1)]   # C3, A2
}

# Possible knight moves (L-shape)
knight_moves = [
    (2, 1), (2, -1), (-2, 1), (-2, -1),
    (1, 2), (1, -2), (-1, 2), (-1, -2)
]

# Function to check if a position is within the board
def is_within_board(x, y):
    return 0 <= x < 3 and 0 <= y < 3

# Function to perform BFS to find the sequence of moves
def bfs_knight_swap():
    # Queue for BFS: (current_positions, move_sequence, turn)
    queue = deque([(initial_positions, [], 'w')])
    visited = set()

    while queue:
        current_positions, move_sequence, turn = queue.popleft()
        current_key = (tuple(current_positions['w']), tuple(current_positions['B']), turn)

        # Check if we reached the target positions
        if current_positions == target_positions:
            return move_sequence

        if current_key in visited:
            continue
        visited.add(current_key)

        # Get the current knight positions for the turn
        knights = current_positions[turn]

        # Try all possible moves for each knight
        for i, (x, y) in enumerate(knights):
            for dx, dy in knight_moves:
                nx, ny = x + dx, y + dy
                if is_within_board(nx, ny) and (nx, ny) not in knights and (nx, ny) not in current_positions['w' if turn == 'B' else 'B']:
                    # Create new positions after the move
                    new_positions = {k: list(v) for k, v in current_positions.items()}
                    new_positions[turn][i] = (nx, ny)
                    new_move_sequence = move_sequence + [f"{turn},{chr(65 + x)}{3 - y},{chr(65 + nx)}{3 - ny}"]
                    # Alternate turn
                    next_turn = 'B' if turn == 'w' else 'w'
                    queue.append((new_positions, new_move_sequence, next_turn))

    return "No"

# Find the solution
solution = bfs_knight_swap()
print(solution)