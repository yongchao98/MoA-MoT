from collections import deque

# Define the board and initial positions
initial_positions = {
    'w': [(0, 0), (1, 0)],  # A1, B1
    'B': [(2, 2), (3, 1)]   # C3, D2
}

# Target positions
target_positions = {
    'w': [(2, 2), (3, 1)],  # C3, D2
    'B': [(0, 0), (1, 0)]   # A1, B1
}

# Possible knight moves
knight_moves = [
    (2, 1), (2, -1), (-2, 1), (-2, -1),
    (1, 2), (1, -2), (-1, 2), (-1, -2)
]

# Convert board positions to a tuple for easy comparison
def board_to_tuple(positions):
    return tuple(sorted(positions['w'] + positions['B']))

# Check if a position is within the board
def is_within_board(x, y):
    return 0 <= x < 4 and 0 <= y < 3

# BFS to find the sequence of moves
def find_moves():
    initial_state = board_to_tuple(initial_positions)
    target_state = board_to_tuple(target_positions)
    
    queue = deque([(initial_positions, [])])
    visited = set()
    visited.add(initial_state)
    
    while queue:
        current_positions, moves = queue.popleft()
        
        # Check if we reached the target state
        if board_to_tuple(current_positions) == target_state:
            return moves
        
        # Determine whose turn it is
        turn = 'w' if len(moves) % 2 == 0 else 'B'
        
        # Try all possible moves for the current player's knights
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
                        new_move = f"{turn},{chr(65+x)}{3-y},{chr(65+nx)}{3-ny}"
                        queue.append((new_positions, moves + [new_move]))
    
    return "No"

# Find and print the moves
result = find_moves()
print(result)