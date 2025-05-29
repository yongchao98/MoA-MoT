from collections import deque

# Define the initial positions of the knights
initial_white_knights = {'C1': 'w', 'D1': 'w'}
initial_black_knights = {'C3': 'B', 'B1': 'B'}
empty_squares = {'A1', 'A2', 'B2', 'D2'}

# Define the target positions for the knights
target_white_positions = {'C3', 'B1'}
target_black_positions = {'C1', 'D1'}

# Define the possible moves for a knight from a given position
def knight_moves(position):
    col, row = position
    moves = []
    col_idx = ord(col) - ord('A')
    row_idx = int(row) - 1
    potential_moves = [
        (col_idx + 2, row_idx + 1), (col_idx + 2, row_idx - 1),
        (col_idx - 2, row_idx + 1), (col_idx - 2, row_idx - 1),
        (col_idx + 1, row_idx + 2), (col_idx + 1, row_idx - 2),
        (col_idx - 1, row_idx + 2), (col_idx - 1, row_idx - 2)
    ]
    for c, r in potential_moves:
        if 0 <= c < 4 and 0 <= r < 3:
            moves.append((chr(c + ord('A')), str(r + 1)))
    return moves

# Check if a move is valid
def is_valid_move(start, end, occupied):
    return end in knight_moves(start) and end not in occupied

# Simulate the moves using BFS
def simulate_moves_bfs():
    # Initial state
    initial_state = (initial_white_knights, initial_black_knights, [])
    queue = deque([initial_state])
    visited = set()

    while queue:
        w_positions, b_positions, move_list = queue.popleft()
        state_key = (frozenset(w_positions.items()), frozenset(b_positions.items()))
        
        if state_key in visited:
            continue
        visited.add(state_key)

        # Check if the target is reached
        if set(w_positions.keys()) == target_white_positions and set(b_positions.keys()) == target_black_positions:
            return move_list

        # Generate new states
        occupied = set(w_positions.keys()).union(b_positions.keys()).union(empty_squares)

        # White's turn
        for w_pos in list(w_positions.keys()):
            for move in knight_moves(w_pos):
                if is_valid_move(w_pos, move, occupied):
                    # Make the move
                    new_w_positions = w_positions.copy()
                    new_w_positions[move] = new_w_positions.pop(w_pos)
                    new_occupied = occupied.copy()
                    new_occupied.remove(w_pos)
                    new_occupied.add(move)
                    new_move_list = move_list + [f"w,{w_pos},{move}"]

                    # Black's turn
                    for b_pos in list(b_positions.keys()):
                        for b_move in knight_moves(b_pos):
                            if is_valid_move(b_pos, b_move, new_occupied):
                                # Make the move
                                new_b_positions = b_positions.copy()
                                new_b_positions[b_move] = new_b_positions.pop(b_pos)
                                new_occupied.remove(b_pos)
                                new_occupied.add(b_move)
                                final_move_list = new_move_list + [f"B,{b_pos},{b_move}"]

                                # Add new state to the queue
                                queue.append((new_w_positions, new_b_positions, final_move_list))

    return None

# Get the result
result = simulate_moves_bfs()
if result:
    print(result)
else:
    print("No")