from collections import deque

# Define the initial and target positions
initial_positions = {'w': [(1, 2), (2, 2)], 'B': [(3, 2), (2, 0)]}
target_positions = {'w': [(3, 2), (2, 0)], 'B': [(1, 2), (2, 2)]}

# Define the board size and possible knight moves
board_size = (4, 3)
knight_moves = [(2, 1), (2, -1), (-2, 1), (-2, -1), (1, 2), (1, -2), (-1, 2), (-1, -2)]

# Convert board positions to a set for quick lookup
def positions_to_set(positions):
    return set(positions['w'] + positions['B'])

# Check if a position is within the board
def is_within_board(pos):
    return 0 <= pos[0] < board_size[0] and 0 <= pos[1] < board_size[1]

# Generate possible moves for a knight
def generate_moves(pos, occupied):
    moves = []
    for move in knight_moves:
        new_pos = (pos[0] + move[0], pos[1] + move[1])
        if is_within_board(new_pos) and new_pos not in occupied:
            moves.append(new_pos)
    return moves

# BFS to find the sequence of moves
def bfs(initial_positions, target_positions):
    queue = deque([(initial_positions, [], 'w')])
    visited = set()
    visited.add(frozenset(positions_to_set(initial_positions)))

    while queue:
        current_positions, moves, current_color = queue.popleft()
        current_set = positions_to_set(current_positions)

        # Check if we reached the target
        if current_set == positions_to_set(target_positions):
            return moves

        # Alternate moves between white and black
        next_color = 'B' if current_color == 'w' else 'w'
        for i, pos in enumerate(current_positions[current_color]):
            # Generate possible moves for the current knight
            for new_pos in generate_moves(pos, current_set):
                new_positions = {k: v[:] for k, v in current_positions.items()}
                new_positions[current_color][i] = new_pos
                new_set = positions_to_set(new_positions)

                if frozenset(new_set) not in visited:
                    visited.add(frozenset(new_set))
                    new_moves = moves + [f"{current_color},{chr(65 + pos[0])}{3 - pos[1]},{chr(65 + new_pos[0])}{3 - new_pos[1]}"]
                    queue.append((new_positions, new_moves, next_color))

    return "No"

# Find the solution
solution = bfs(initial_positions, target_positions)
print(solution)