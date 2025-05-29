from collections import deque

# Define the board and initial positions
initial_positions = {'w': [(3, 4), (2, 3)], 'B': [(1, 2), (1, 3)]}
target_positions = {'w': [(1, 2), (1, 3)], 'B': [(3, 4), (2, 3)]}

# Define knight moves (L-shape)
knight_moves = [(2, 1), (2, -1), (-2, 1), (-2, -1), (1, 2), (1, -2), (-1, 2), (-1, -2)]

# Convert board positions to a tuple for easy comparison
def board_to_tuple(positions):
    return tuple(sorted(positions['w'] + positions['B']))

# Check if a position is within the board
def is_within_board(x, y):
    return 1 <= x <= 3 and 1 <= y <= 4

# BFS to find the sequence of moves
def bfs_knight_swap():
    queue = deque([(initial_positions, [])])
    visited = set()
    visited.add(board_to_tuple(initial_positions))
    
    while queue:
        current_positions, moves = queue.popleft()
        
        # Check if we reached the target configuration
        if board_to_tuple(current_positions) == board_to_tuple(target_positions):
            return moves
        
        # Determine whose turn it is
        current_color = 'w' if len(moves) % 2 == 0 else 'B'
        
        # Try all possible moves for the current color
        for i, (x, y) in enumerate(current_positions[current_color]):
            for dx, dy in knight_moves:
                nx, ny = x + dx, y + dy
                if is_within_board(nx, ny) and (nx, ny) not in current_positions['w'] and (nx, ny) not in current_positions['B']:
                    # Make the move
                    new_positions = {color: list(pos) for color, pos in current_positions.items()}
                    new_positions[current_color][i] = (nx, ny)
                    new_positions_tuple = board_to_tuple(new_positions)
                    
                    if new_positions_tuple not in visited:
                        visited.add(new_positions_tuple)
                        new_move = f"{current_color},{chr(64+y)}{x},{chr(64+ny)}{nx}"
                        queue.append((new_positions, moves + [new_move]))
    
    return "No"

# Run the BFS to find the solution
result = bfs_knight_swap()
print(result)