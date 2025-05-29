from collections import deque

# Define the board and initial positions
initial_positions = {
    'w': [(1, 1), (0, 1)],  # White knights at B2, B1
    'B': [(2, 3), (0, 0)]   # Black knights at D3, A1
}

# Target positions
target_positions = {
    'w': [(2, 3), (0, 0)],  # White knights target D3, A1
    'B': [(1, 1), (0, 1)]   # Black knights target B2, B1
}

# Knight moves (L-shape)
knight_moves = [
    (2, 1), (2, -1), (-2, 1), (-2, -1),
    (1, 2), (1, -2), (-1, 2), (-1, -2)
]

# Convert board positions to a tuple for easy comparison
def board_to_tuple(positions):
    return tuple(sorted(positions['w'] + positions['B']))

# Check if a position is within the board
def is_within_board(x, y):
    return 0 <= x < 3 and 0 <= y < 4

# BFS to find the sequence of moves
def bfs_knight_swap():
    # Initial state
    start_state = board_to_tuple(initial_positions)
    target_state = board_to_tuple(target_positions)
    
    # Queue for BFS
    queue = deque([(initial_positions, [])])
    visited = set()
    visited.add(start_state)
    
    while queue:
        current_positions, moves = queue.popleft()
        
        # Check if we reached the target state
        if board_to_tuple(current_positions) == target_state:
            return moves
        
        # Determine whose turn it is
        turn = 'B' if len(moves) % 2 == 0 else 'w'
        
        # Try all possible moves for the current player
        for i, (x, y) in enumerate(current_positions[turn]):
            for dx, dy in knight_moves:
                nx, ny = x + dx, y + dy
                if is_within_board(nx, ny) and (nx, ny) not in current_positions['w'] and (nx, ny) not in current_positions['B']:
                    # Make the move
                    new_positions = { 'w': current_positions['w'][:], 'B': current_positions['B'][:] }
                    new_positions[turn][i] = (nx, ny)
                    new_state = board_to_tuple(new_positions)
                    
                    if new_state not in visited:
                        visited.add(new_state)
                        new_move = f"{turn},{chr(65+y)}{3-x},{chr(65+ny)}{3-nx}"
                        queue.append((new_positions, moves + [new_move]))
    
    return "No"

# Execute the BFS to find the solution
result = bfs_knight_swap()
print(result)