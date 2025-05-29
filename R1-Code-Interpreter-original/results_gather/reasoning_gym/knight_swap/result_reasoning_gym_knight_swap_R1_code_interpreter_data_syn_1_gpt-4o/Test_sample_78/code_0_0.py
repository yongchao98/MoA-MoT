from collections import deque

# Define the initial and target positions
initial_positions = {'w': [(2, 2), (1, 2)], 'B': [(0, 0), (1, 1)]}
target_positions = {'w': [(0, 0), (1, 1)], 'B': [(2, 2), (1, 2)]}

# Define the knight moves
knight_moves = [(2, 1), (2, -1), (-2, 1), (-2, -1), (1, 2), (1, -2), (-1, 2), (-1, -2)]

# Convert board positions to a tuple for easy comparison
def board_to_tuple(positions):
    return tuple(sorted(positions['w'] + positions['B']))

# Check if the current positions match the target positions
def is_target(positions):
    return sorted(positions['w']) == sorted(target_positions['w']) and sorted(positions['B']) == sorted(target_positions['B'])

# BFS to find the sequence of moves
def find_moves():
    queue = deque([(initial_positions, [])])
    visited = set()
    visited.add(board_to_tuple(initial_positions))
    
    while queue:
        current_positions, moves = queue.popleft()
        
        if is_target(current_positions):
            return moves
        
        # Determine whose turn it is
        current_color = 'w' if len(moves) % 2 == 0 else 'B'
        
        for i, (x, y) in enumerate(current_positions[current_color]):
            for dx, dy in knight_moves:
                nx, ny = x + dx, y + dy
                if 0 <= nx < 4 and 0 <= ny < 3 and (nx, ny) not in current_positions['w'] and (nx, ny) not in current_positions['B']:
                    new_positions = {'w': current_positions['w'][:], 'B': current_positions['B'][:]}
                    new_positions[current_color][i] = (nx, ny)
                    new_positions_tuple = board_to_tuple(new_positions)
                    
                    if new_positions_tuple not in visited:
                        visited.add(new_positions_tuple)
                        new_move = f"{current_color},{chr(65 + x)}{3 - y},{chr(65 + nx)}{3 - ny}"
                        queue.append((new_positions, moves + [new_move]))
    
    return "No"

# Find and print the moves
result = find_moves()
print(result)