from collections import deque

# Define the board and initial positions
initial_positions = {
    'w': [(1, 2), (3, 2)],  # B3, D3
    'B': [(0, 1), (2, 0)]   # A2, C1
}

# Target positions for the swap
target_positions = {
    'w': [(0, 1), (2, 0)],  # A2, C1
    'B': [(1, 2), (3, 2)]   # B3, D3
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
            board[y][x] = color
    return ''.join(''.join(row) for row in board)

# Check if a position is within the board
def is_within_board(x, y):
    return 0 <= x < 4 and 0 <= y < 3

# Perform BFS to find the sequence of moves
def find_moves():
    initial_state = board_to_string(initial_positions)
    target_state = board_to_string(target_positions)
    
    queue = deque([(initial_positions, [], 'w')])  # (current positions, moves, next color)
    visited = set()
    visited.add(initial_state)
    
    while queue:
        current_positions, moves, next_color = queue.popleft()
        
        # Check if we reached the target state
        if board_to_string(current_positions) == target_state:
            return moves
        
        # Get the positions of the knights of the current color
        for i, (x, y) in enumerate(current_positions[next_color]):
            for dx, dy in knight_moves:
                nx, ny = x + dx, y + dy
                if is_within_board(nx, ny) and (nx, ny) not in current_positions['w'] and (nx, ny) not in current_positions['B']:
                    # Create a new state with the knight moved
                    new_positions = {color: list(pos_list) for color, pos_list in current_positions.items()}
                    new_positions[next_color][i] = (nx, ny)
                    
                    # Convert to string to check if visited
                    new_state = board_to_string(new_positions)
                    if new_state not in visited:
                        visited.add(new_state)
                        new_moves = moves + [f"{next_color},{chr(x + 65)}{3 - y},{chr(nx + 65)}{3 - ny}"]
                        queue.append((new_positions, new_moves, 'B' if next_color == 'w' else 'w'))
    
    return "No"

# Find the sequence of moves
result = find_moves()
print(result)