from collections import deque

# Define the board and initial positions
initial_positions = {
    'w': [(0, 1), (0, 2)],  # B3, C3
    'B': [(0, 0), (1, 3)]   # A3, D2
}

# Target positions
target_positions = {
    'w': [(0, 0), (1, 3)],  # A3, D2
    'B': [(0, 1), (0, 2)]   # B3, C3
}

# Possible knight moves (L-shape)
knight_moves = [
    (2, 1), (2, -1), (-2, 1), (-2, -1),
    (1, 2), (1, -2), (-1, 2), (-1, -2)
]

# Convert board positions to a string for easy comparison
def board_to_string(positions):
    board = [['.' for _ in range(4)] for _ in range(3)]
    for color, pos_list in positions.items():
        for x, y in pos_list:
            board[x][y] = color
    return ''.join(''.join(row) for row in board)

# Check if a position is within the board
def is_within_board(x, y):
    return 0 <= x < 3 and 0 <= y < 4

# BFS to find the sequence of moves
def find_moves():
    initial_state = board_to_string(initial_positions)
    target_state = board_to_string(target_positions)
    queue = deque([(initial_positions, [], 'w')])  # (positions, moves, current_turn)
    visited = set()
    visited.add(initial_state)

    while queue:
        current_positions, moves, turn = queue.popleft()
        current_state = board_to_string(current_positions)

        if current_state == target_state:
            return moves

        # Get the positions of the current player's knights
        for i, (x, y) in enumerate(current_positions[turn]):
            for dx, dy in knight_moves:
                nx, ny = x + dx, y + dy
                if is_within_board(nx, ny) and (nx, ny) not in current_positions['w'] and (nx, ny) not in current_positions['B']:
                    # Make the move
                    new_positions = {color: list(pos_list) for color, pos_list in current_positions.items()}
                    new_positions[turn][i] = (nx, ny)
                    new_state = board_to_string(new_positions)
                    if new_state not in visited:
                        visited.add(new_state)
                        new_moves = moves + [f"{turn},{chr(65+y)}{3-x},{chr(65+ny)}{3-nx}"]
                        queue.append((new_positions, new_moves, 'B' if turn == 'w' else 'w'))

    return "No"

# Find the moves
result = find_moves()
print(result)