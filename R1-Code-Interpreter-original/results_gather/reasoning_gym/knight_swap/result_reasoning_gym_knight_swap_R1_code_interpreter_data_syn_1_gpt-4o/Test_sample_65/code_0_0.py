from collections import deque

# Define the board and initial positions
initial_positions = {
    'w': [(0, 2), (1, 0)],  # A3, B1
    'B': [(2, 2), (3, 1)]   # C3, D2
}

# Target positions
target_positions = {
    'w': [(2, 2), (3, 1)],  # C3, D2
    'B': [(0, 2), (1, 0)]   # A3, B1
}

# Possible knight moves (L-shape)
knight_moves = [
    (2, 1), (2, -1), (-2, 1), (-2, -1),
    (1, 2), (1, -2), (-1, 2), (-1, -2)
]

# Convert board coordinates to chess notation
def to_chess_notation(x, y):
    return chr(ord('A') + x) + str(y + 1)

# Check if a position is within the board
def is_within_board(x, y):
    return 0 <= x < 4 and 0 <= y < 3

# BFS to find the sequence of moves
def bfs_knight_swap():
    queue = deque([(initial_positions, [])])
    visited = set()
    visited.add((tuple(initial_positions['w']), tuple(initial_positions['B'])))
    
    while queue:
        positions, moves = queue.popleft()
        
        # Check if we reached the target positions
        if positions == target_positions:
            return moves
        
        # Determine whose turn it is
        current_color = 'B' if len(moves) % 2 == 0 else 'w'
        
        # Try all possible moves for the current color
        for i, (x, y) in enumerate(positions[current_color]):
            for dx, dy in knight_moves:
                nx, ny = x + dx, y + dy
                if is_within_board(nx, ny) and (nx, ny) not in positions['w'] and (nx, ny) not in positions['B']:
                    # Make the move
                    new_positions = {color: list(pos) for color, pos in positions.items()}
                    new_positions[current_color][i] = (nx, ny)
                    new_positions_tuple = (tuple(new_positions['w']), tuple(new_positions['B']))
                    
                    if new_positions_tuple not in visited:
                        visited.add(new_positions_tuple)
                        move_notation = f"{current_color},{to_chess_notation(x, y)},{to_chess_notation(nx, ny)}"
                        queue.append((new_positions, moves + [move_notation]))
    
    return "No"

# Execute the BFS to find the solution
result = bfs_knight_swap()
print(result)