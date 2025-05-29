from collections import deque

# Define the initial board state
initial_state = {
    'w': [(0, 1), (3, 1)],  # A2, D2
    'B': [(0, 2), (1, 2)]   # A3, B3
}

# Define the target board state
target_state = {
    'w': [(0, 2), (1, 2)],  # A3, B3
    'B': [(0, 1), (3, 1)]   # A2, D2
}

# Define the board size
board_size = (4, 3)  # 4 columns (A-D), 3 rows (1-3)

# Define knight moves (L-shape)
knight_moves = [
    (2, 1), (2, -1), (-2, 1), (-2, -1),
    (1, 2), (1, -2), (-1, 2), (-1, -2)
]

# Convert position to board coordinates
def pos_to_coord(pos):
    col, row = pos
    return ord(col) - ord('A'), int(row) - 1

# Convert board coordinates to position
def coord_to_pos(coord):
    x, y = coord
    return chr(x + ord('A')) + str(y + 1)

# Check if a position is within the board
def is_within_board(x, y):
    return 0 <= x < board_size[0] and 0 <= y < board_size[1]

# Generate possible moves for a knight
def generate_moves(knight_pos, occupied_positions):
    moves = []
    for dx, dy in knight_moves:
        new_x, new_y = knight_pos[0] + dx, knight_pos[1] + dy
        if is_within_board(new_x, new_y) and (new_x, new_y) not in occupied_positions:
            moves.append((new_x, new_y))
    return moves

# BFS to find the sequence of moves
def bfs_knight_swap(initial_state, target_state):
    queue = deque([(initial_state, [], 'w')])  # (current_state, moves, current_turn)
    visited = set()

    while queue:
        current_state, moves, current_turn = queue.popleft()
        state_key = (tuple(current_state['w']), tuple(current_state['B']), current_turn)
        
        if state_key in visited:
            continue
        visited.add(state_key)

        # Check if the current state matches the target state
        if current_state == target_state:
            return moves

        # Get the positions of the current player's knights
        current_knights = current_state[current_turn]
        other_knights = current_state['B' if current_turn == 'w' else 'w']
        occupied_positions = set(current_knights + other_knights)

        # Generate all possible moves for the current player's knights
        for i, knight_pos in enumerate(current_knights):
            possible_moves = generate_moves(knight_pos, occupied_positions)
            for new_pos in possible_moves:
                # Create a new state with the knight moved
                new_knights = current_knights[:i] + [new_pos] + current_knights[i+1:]
                new_state = {
                    current_turn: new_knights,
                    'B' if current_turn == 'w' else 'w': other_knights
                }
                # Record the move
                move = f"{current_turn},{coord_to_pos(knight_pos)},{coord_to_pos(new_pos)}"
                # Add the new state to the queue
                queue.append((new_state, moves + [move], 'B' if current_turn == 'w' else 'w'))

    return "No"

# Find the solution
solution = bfs_knight_swap(initial_state, target_state)
print(solution)