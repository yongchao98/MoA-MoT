from collections import deque

# Define the board and initial positions
initial_positions = {'w': [(1, 2), (0, 1)], 'B': [(2, 3), (1, 0)]}
target_positions = {'w': [(2, 3), (1, 0)], 'B': [(1, 2), (0, 1)]}

# Define the knight moves
knight_moves = [(2, 1), (2, -1), (-2, 1), (-2, -1), (1, 2), (1, -2), (-1, 2), (-1, -2)]

# Convert board positions to a tuple for easy comparison
def board_to_tuple(positions):
    return tuple(sorted(positions['w'] + positions['B']))

# Check if a position is within the board
def is_within_board(x, y):
    return 0 <= x < 3 and 0 <= y < 4

# BFS to find the sequence of moves
def bfs_knight_swap():
    queue = deque([(initial_positions, [])])
    visited = set()
    visited.add(board_to_tuple(initial_positions))
    
    while queue:
        current_positions, moves = queue.popleft()
        
        # Check if we reached the target positions
        if board_to_tuple(current_positions) == board_to_tuple(target_positions):
            return moves
        
        # Determine whose turn it is
        turn = 'w' if len(moves) % 2 == 0 else 'B'
        
        # Try all possible moves for the current player's knights
        for i, (x, y) in enumerate(current_positions[turn]):
            for dx, dy in knight_moves:
                nx, ny = x + dx, y + dy
                if is_within_board(nx, ny) and (nx, ny) not in current_positions['w'] and (nx, ny) not in current_positions['B']:
                    # Make the move
                    new_positions = {k: v[:] for k, v in current_positions.items()}
                    new_positions[turn][i] = (nx, ny)
                    new_positions_tuple = board_to_tuple(new_positions)
                    
                    if new_positions_tuple not in visited:
                        visited.add(new_positions_tuple)
                        move_str = f"{turn},{chr(65 + y)}{3 - x},{chr(65 + ny)}{3 - nx}"
                        queue.append((new_positions, moves + [move_str]))
    
    return "No"

# Execute the BFS to find the solution
result = bfs_knight_swap()
print(result)