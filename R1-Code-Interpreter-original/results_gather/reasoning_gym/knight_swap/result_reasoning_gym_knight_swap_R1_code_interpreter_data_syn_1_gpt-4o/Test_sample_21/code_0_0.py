from collections import deque

# Define the initial and target positions
initial_positions = {'w': [(2, 1), (3, 2)], 'B': [(1, 3), (3, 3)]}
target_positions = {'w': [(1, 3), (3, 3)], 'B': [(2, 1), (3, 2)]}

# Define the board size
board_size = (4, 3)

# Define knight moves (L-shape)
knight_moves = [(2, 1), (2, -1), (-2, 1), (-2, -1),
                (1, 2), (1, -2), (-1, 2), (-1, -2)]

# Function to check if a position is within the board
def is_within_board(x, y):
    return 1 <= x <= board_size[0] and 1 <= y <= board_size[1]

# Function to generate possible moves for a knight
def generate_moves(position, occupied_positions):
    x, y = position
    possible_moves = []
    for dx, dy in knight_moves:
        new_x, new_y = x + dx, y + dy
        if is_within_board(new_x, new_y) and (new_x, new_y) not in occupied_positions:
            possible_moves.append((new_x, new_y))
    return possible_moves

# BFS to find the sequence of moves
def bfs_knight_swap():
    # Initial state
    initial_state = (tuple(initial_positions['w']), tuple(initial_positions['B']), 'w')
    queue = deque([(initial_state, [])])
    visited = set([initial_state])

    while queue:
        (white_positions, black_positions, turn), moves = queue.popleft()

        # Check if we reached the target configuration
        if (white_positions == tuple(target_positions['w']) and
            black_positions == tuple(target_positions['B'])):
            return moves

        # Determine the current player's knights and the other player's knights
        current_knights = white_positions if turn == 'w' else black_positions
        other_knights = black_positions if turn == 'w' else white_positions

        # Generate all possible moves for the current player's knights
        for i, knight_pos in enumerate(current_knights):
            occupied_positions = set(white_positions + black_positions)
            possible_moves = generate_moves(knight_pos, occupied_positions)

            for move in possible_moves:
                new_knights = list(current_knights)
                new_knights[i] = move
                new_knights = tuple(new_knights)

                # Create the new state
                if turn == 'w':
                    new_state = (new_knights, black_positions, 'B')
                else:
                    new_state = (white_positions, new_knights, 'w')

                if new_state not in visited:
                    visited.add(new_state)
                    new_move = f"{turn},{chr(65 + knight_pos[0] - 1)}{knight_pos[1]},{chr(65 + move[0] - 1)}{move[1]}"
                    queue.append((new_state, moves + [new_move]))

    return "No"

# Execute the BFS to find the solution
result = bfs_knight_swap()
print(result)