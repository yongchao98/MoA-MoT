from collections import deque

# Define the board and initial positions
initial_positions = {'w': [(1, 2), (2, 0)], 'B': [(0, 1), (2, 2)]}
target_positions = {'w': [(0, 1), (2, 2)], 'B': [(1, 2), (2, 0)]}
board_size = (3, 4)

# Define knight moves
knight_moves = [(2, 1), (2, -1), (-2, 1), (-2, -1), (1, 2), (1, -2), (-1, 2), (-1, -2)]

# Convert position to board notation
def pos_to_notation(pos):
    return chr(pos[1] + ord('A')) + str(pos[0] + 1)

# Check if a position is within the board
def is_within_board(pos):
    return 0 <= pos[0] < board_size[0] and 0 <= pos[1] < board_size[1]

# BFS to find the sequence of moves
def find_moves():
    queue = deque([(initial_positions, [])])
    visited = set()
    visited.add(tuple(sorted(initial_positions['w'] + initial_positions['B'])))

    while queue:
        current_positions, moves = queue.popleft()

        # Check if we reached the target positions
        if sorted(current_positions['w']) == sorted(target_positions['w']) and \
           sorted(current_positions['B']) == sorted(target_positions['B']):
            return moves

        # Determine whose turn it is
        turn = 'w' if len(moves) % 2 == 0 else 'B'

        # Try all possible moves for the current player
        for i, knight_pos in enumerate(current_positions[turn]):
            for move in knight_moves:
                new_pos = (knight_pos[0] + move[0], knight_pos[1] + move[1])
                if is_within_board(new_pos) and new_pos not in current_positions['w'] and new_pos not in current_positions['B']:
                    new_positions = {k: list(v) for k, v in current_positions.items()}
                    new_positions[turn][i] = new_pos
                    new_positions_tuple = tuple(sorted(new_positions['w'] + new_positions['B']))
                    if new_positions_tuple not in visited:
                        visited.add(new_positions_tuple)
                        new_move = f"{turn},{pos_to_notation(knight_pos)},{pos_to_notation(new_pos)}"
                        queue.append((new_positions, moves + [new_move]))

    return "No"

# Find and print the moves
result = find_moves()
print(result)