from collections import deque

# Define the board and initial positions
initial_positions = {
    'w': [(1, 2), (2, 0)],  # B3, C1
    'B': [(0, 1), (2, 2)]   # A2, C3
}

# Target positions
target_positions = {
    'w': [(0, 1), (2, 2)],  # A2, C3
    'B': [(1, 2), (2, 0)]   # B3, C1
}

# Possible knight moves
knight_moves = [
    (2, 1), (2, -1), (-2, 1), (-2, -1),
    (1, 2), (1, -2), (-1, 2), (-1, -2)
]

# Convert position to board notation
def pos_to_notation(pos):
    return chr(pos[0] + ord('A')) + str(pos[1] + 1)

# Check if a position is within the board
def is_within_board(pos):
    return 0 <= pos[0] < 4 and 0 <= pos[1] < 3

# BFS to find the sequence of moves
def bfs_knight_swap():
    queue = deque([(initial_positions, [])])
    visited = set()
    visited.add((tuple(initial_positions['w']), tuple(initial_positions['B'])))

    while queue:
        current_positions, moves = queue.popleft()

        # Check if we reached the target positions
        if current_positions == target_positions:
            return moves

        # Determine whose turn it is
        turn = 'w' if len(moves) % 2 == 0 else 'B'

        # Try all possible moves for the current player
        for i, knight_pos in enumerate(current_positions[turn]):
            for move in knight_moves:
                new_pos = (knight_pos[0] + move[0], knight_pos[1] + move[1])

                if is_within_board(new_pos) and new_pos not in current_positions['w'] and new_pos not in current_positions['B']:
                    # Create new positions
                    new_positions = {
                        'w': list(current_positions['w']),
                        'B': list(current_positions['B'])
                    }
                    new_positions[turn][i] = new_pos

                    # Convert to tuple for immutability
                    new_positions_tuple = (tuple(new_positions['w']), tuple(new_positions['B']))

                    if new_positions_tuple not in visited:
                        visited.add(new_positions_tuple)
                        new_move = f"{turn},{pos_to_notation(knight_pos)},{pos_to_notation(new_pos)}"
                        queue.append((new_positions, moves + [new_move]))

    return "No"

# Execute the BFS to find the solution
result = bfs_knight_swap()
print(result)