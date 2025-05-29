from collections import deque

# Define the initial and target positions
initial_positions = {'w': [(3, 4), (2, 2)], 'B': [(3, 2), (3, 3)]}
target_positions = {'w': [(3, 2), (3, 3)], 'B': [(3, 4), (2, 2)]}

# Define the board size and knight moves
board_size = (3, 4)
knight_moves = [(2, 1), (2, -1), (-2, 1), (-2, -1), (1, 2), (1, -2), (-1, 2), (-1, -2)]

# Convert board positions to a more manageable format
def pos_to_str(pos):
    return f"{chr(64 + pos[1])}{pos[0]}"

def is_within_board(pos):
    return 1 <= pos[0] <= board_size[0] and 1 <= pos[1] <= board_size[1]

# BFS to find the sequence of moves
def find_moves():
    queue = deque([(initial_positions, [])])
    visited = set()

    while queue:
        current_positions, moves = queue.popleft()
        current_state = (tuple(current_positions['w']), tuple(current_positions['B']))

        if current_state in visited:
            continue
        visited.add(current_state)

        # Check if the target is reached
        if current_positions == target_positions:
            return moves

        # Generate possible moves for white and black knights
        for color in ['w', 'B']:
            for i, knight_pos in enumerate(current_positions[color]):
                for move in knight_moves:
                    new_pos = (knight_pos[0] + move[0], knight_pos[1] + move[1])
                    if is_within_board(new_pos) and new_pos not in current_positions['w'] and new_pos not in current_positions['B']:
                        new_positions = {k: v[:] for k, v in current_positions.items()}
                        new_positions[color][i] = new_pos
                        new_moves = moves + [f"{color},{pos_to_str(knight_pos)},{pos_to_str(new_pos)}"]
                        queue.append((new_positions, new_moves))

    return "No"

# Find and print the sequence of moves
result = find_moves()
print(result)