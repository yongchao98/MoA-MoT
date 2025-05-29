from collections import deque

# Define the board and initial positions
initial_positions = {
    'w': [(3, 3), (1, 2)],  # D3, C1
    'B': [(3, 1), (2, 2)]   # B3, C2
}

# Target positions
target_positions = {
    'w': [(3, 1), (2, 2)],  # B3, C2
    'B': [(3, 3), (1, 2)]   # D3, C1
}

# Possible knight moves (L-shape)
knight_moves = [
    (2, 1), (2, -1), (-2, 1), (-2, -1),
    (1, 2), (1, -2), (-1, 2), (-1, -2)
]

# Convert board positions to a tuple for easy comparison
def board_to_tuple(positions):
    return tuple(sorted(positions['w'] + positions['B']))

# Check if the current positions match the target positions
def is_target(positions):
    return sorted(positions['w']) == sorted(target_positions['w']) and \
           sorted(positions['B']) == sorted(target_positions['B'])

# BFS to find the sequence of moves
def find_moves():
    # Queue for BFS: (current positions, move sequence, current player)
    queue = deque([(initial_positions, [], 'w')])
    visited = set()
    visited.add(board_to_tuple(initial_positions))

    while queue:
        current_positions, moves, current_player = queue.popleft()

        # Check if we reached the target
        if is_target(current_positions):
            return moves

        # Get the positions of the current player's knights
        for i, (x, y) in enumerate(current_positions[current_player]):
            # Try all possible knight moves
            for dx, dy in knight_moves:
                nx, ny = x + dx, y + dy
                # Check if the move is within bounds and to an empty square
                if 1 <= nx <= 3 and 1 <= ny <= 4 and (nx, ny) not in current_positions['w'] and (nx, ny) not in current_positions['B']:
                    # Create new positions after the move
                    new_positions = {
                        'w': current_positions['w'][:],
                        'B': current_positions['B'][:]
                    }
                    new_positions[current_player][i] = (nx, ny)
                    new_positions_tuple = board_to_tuple(new_positions)

                    # If not visited, add to queue
                    if new_positions_tuple not in visited:
                        visited.add(new_positions_tuple)
                        # Record the move
                        move = f"{current_player},{chr(64 + y)}{x},{chr(64 + ny)}{nx}"
                        queue.append((new_positions, moves + [move], 'B' if current_player == 'w' else 'w'))

    return "No"

# Find the moves
result = find_moves()
print(result)