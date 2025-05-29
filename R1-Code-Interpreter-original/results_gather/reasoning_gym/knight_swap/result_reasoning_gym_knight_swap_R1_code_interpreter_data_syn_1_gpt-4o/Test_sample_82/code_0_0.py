from collections import deque

# Define the board and initial positions
board = {
    'A1': 'B', 'B1': 'B', 'C1': '.', 'D1': '.',
    'A2': 'w', 'B2': '.', 'C2': '.', 'D2': '.',
    'A3': '.', 'B3': 'w', 'C3': '.', 'D3': '.'
}

# Define target positions
target_positions = {
    'w': ['A1', 'B1'],
    'B': ['B3', 'A2']
}

# Define knight moves (L-shape)
knight_moves = [
    (2, 1), (2, -1), (-2, 1), (-2, -1),
    (1, 2), (1, -2), (-1, 2), (-1, -2)
]

# Convert board position to coordinates
def pos_to_coord(pos):
    col, row = pos
    return ord(col) - ord('A'), int(row) - 1

# Convert coordinates to board position
def coord_to_pos(coord):
    col, row = coord
    return chr(col + ord('A')) + str(row + 1)

# Check if a position is within the board
def is_within_board(coord):
    col, row = coord
    return 0 <= col < 4 and 0 <= row < 3

# BFS to find the sequence of moves
def bfs_knight_swap():
    queue = deque()
    visited = set()
    initial_state = (('B3', 'A2'), ('A1', 'B1'))  # (white_positions, black_positions)
    queue.append((initial_state, [], 'w'))  # (state, moves, current_turn)
    visited.add(initial_state)

    while queue:
        (white_positions, black_positions), moves, current_turn = queue.popleft()

        # Check if the current state is the target state
        if set(white_positions) == set(target_positions['w']) and set(black_positions) == set(target_positions['B']):
            return moves

        # Determine current positions and target positions
        current_positions = white_positions if current_turn == 'w' else black_positions
        other_positions = black_positions if current_turn == 'w' else white_positions

        # Try all possible moves for the current knight
        for i, pos in enumerate(current_positions):
            col, row = pos_to_coord(pos)
            for move in knight_moves:
                new_coord = (col + move[0], row + move[1])
                if is_within_board(new_coord):
                    new_pos = coord_to_pos(new_coord)
                    if new_pos not in current_positions and new_pos not in other_positions and board[new_pos] == '.':
                        # Create new state
                        new_positions = list(current_positions)
                        new_positions[i] = new_pos
                        new_positions = tuple(new_positions)
                        new_state = (new_positions, other_positions) if current_turn == 'w' else (other_positions, new_positions)
                        if new_state not in visited:
                            visited.add(new_state)
                            new_moves = moves + [f"{current_turn},{pos},{new_pos}"]
                            queue.append((new_state, new_moves, 'B' if current_turn == 'w' else 'w'))

    return "No"

# Execute the BFS to find the solution
result = bfs_knight_swap()
print(result)