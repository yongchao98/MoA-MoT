from collections import deque

# Define the board and initial positions
initial_positions = {
    'w': [(1, 1), (3, 2)],  # B1, D2
    'B': [(0, 2), (1, 2)]   # A3, B3
}

# Target positions for each color
target_positions = {
    'w': [(0, 2), (1, 2)],  # A3, B3
    'B': [(1, 1), (3, 2)]   # B1, D2
}

# Possible knight moves
knight_moves = [(2, 1), (1, 2), (-1, 2), (-2, 1), (-2, -1), (-1, -2), (1, -2), (2, -1)]

# Convert board positions to a more manageable form
def pos_to_str(pos):
    return chr(pos[0] + ord('A')) + str(pos[1] + 1)

# Check if a position is within the board
def is_within_board(pos):
    return 0 <= pos[0] < 4 and 0 <= pos[1] < 3

# Convert positions dictionary to a hashable form
def positions_to_hashable(positions):
    return (tuple(sorted(positions['w'])), tuple(sorted(positions['B'])))

# BFS to find the sequence of moves
def bfs_knight_swap():
    queue = deque([(initial_positions, [])])
    visited = set()

    while queue:
        positions, moves = queue.popleft()
        if positions == target_positions:
            return moves

        # Alternate moves between white and black
        current_color = 'w' if len(moves) % 2 == 0 else 'B'
        for i, knight_pos in enumerate(positions[current_color]):
            for move in knight_moves:
                new_pos = (knight_pos[0] + move[0], knight_pos[1] + move[1])
                if is_within_board(new_pos) and new_pos not in positions['w'] and new_pos not in positions['B']:
                    new_positions = {color: list(pos) for color, pos in positions.items()}
                    new_positions[current_color][i] = new_pos
                    new_positions = {color: tuple(pos) for color, pos in new_positions.items()}
                    hashable_positions = positions_to_hashable(new_positions)
                    if hashable_positions not in visited:
                        visited.add(hashable_positions)
                        new_moves = moves + [f"{current_color},{pos_to_str(knight_pos)},{pos_to_str(new_pos)}"]
                        queue.append((new_positions, new_moves))
    return "No"

# Execute the BFS to find the solution
result = bfs_knight_swap()
print(result)