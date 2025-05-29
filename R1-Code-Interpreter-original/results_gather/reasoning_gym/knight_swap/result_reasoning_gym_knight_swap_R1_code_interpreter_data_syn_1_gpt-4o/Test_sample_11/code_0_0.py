from collections import deque

# Define the board and initial positions
initial_positions = {'w': [(1, 2), (0, 1)], 'B': [(3, 1), (2, 0)]}
target_positions = {'w': [(3, 1), (2, 0)], 'B': [(1, 2), (0, 1)]}

# Define knight moves
knight_moves = [(2, 1), (1, 2), (-1, 2), (-2, 1), (-2, -1), (-1, -2), (1, -2), (2, -1)]

# Convert position to board notation
def pos_to_notation(pos):
    return chr(pos[0] + ord('A')) + str(pos[1] + 1)

# Check if a position is within the board
def is_within_board(pos):
    return 0 <= pos[0] < 4 and 0 <= pos[1] < 3

# BFS to find the sequence of moves
def bfs_knight_swap():
    queue = deque([(initial_positions, [], 'w')])  # (positions, moves, current_turn)
    visited = set()

    while queue:
        positions, moves, current_turn = queue.popleft()
        if positions == target_positions:
            return moves

        # Generate possible moves for the current turn
        for i, knight_pos in enumerate(positions[current_turn]):
            for move in knight_moves:
                new_pos = (knight_pos[0] + move[0], knight_pos[1] + move[1])
                if is_within_board(new_pos) and new_pos not in positions['w'] and new_pos not in positions['B']:
                    new_positions = {color: list(pos) for color, pos in positions.items()}
                    new_positions[current_turn][i] = new_pos
                    new_positions = {color: tuple(pos) for color, pos in new_positions.items()}
                    if new_positions not in visited:
                        visited.add(new_positions)
                        new_moves = moves + [f"{current_turn},{pos_to_notation(knight_pos)},{pos_to_notation(new_pos)}"]
                        next_turn = 'B' if current_turn == 'w' else 'w'
                        queue.append((new_positions, new_moves, next_turn))

    return "No"

# Find the solution
solution = bfs_knight_swap()
print(solution)