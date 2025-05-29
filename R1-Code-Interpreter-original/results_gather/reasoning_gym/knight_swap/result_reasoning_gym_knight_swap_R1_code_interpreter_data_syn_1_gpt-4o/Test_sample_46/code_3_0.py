from collections import deque

# Define the initial positions of the knights
initial_white = {'B3': 'w', 'A1': 'w'}
initial_black = {'D2': 'B', 'C1': 'B'}
initial_empty = {'A2', 'D3', 'B1'}

# Define the target positions for the swap
target_white = {'D2', 'C1'}
target_black = {'B3', 'A1'}

# Define the possible moves for a knight
knight_moves = [(2, 1), (1, 2), (-1, 2), (-2, 1), (-2, -1), (-1, -2), (1, -2), (2, -1)]

# Convert board positions to coordinates
def pos_to_coord(pos):
    col, row = pos
    return ord(col) - ord('A'), int(row) - 1

# Convert coordinates back to board positions
def coord_to_pos(coord):
    col, row = coord
    return chr(col + ord('A')) + str(row + 1)

# Check if a position is within the board
def is_within_board(coord):
    col, row = coord
    return 0 <= col < 4 and 0 <= row < 3

# Generate possible moves for a knight
def generate_moves(position, empty_squares):
    start_coord = pos_to_coord(position)
    possible_moves = []
    for move in knight_moves:
        end_coord = (start_coord[0] + move[0], start_coord[1] + move[1])
        if is_within_board(end_coord):
            end_pos = coord_to_pos(end_coord)
            if end_pos in empty_squares:
                possible_moves.append(end_pos)
    return possible_moves

# Simulate the moves using BFS
def simulate_moves():
    # Initial state
    initial_state = (frozenset(initial_white.items()), frozenset(initial_black.items()), frozenset(initial_empty))
    target_state = (frozenset(target_white), frozenset(target_black))

    # Queue for BFS
    queue = deque([(initial_state, [])])
    visited = set()

    while queue:
        (current_white, current_black, current_empty), path = queue.popleft()

        # Check if target state is reached
        if (current_white == target_black) and (current_black == target_white):
            return path

        # Convert sets back to dictionaries for manipulation
        current_white = dict(current_white)
        current_black = dict(current_black)
        current_empty = set(current_empty)

        # Generate moves for black knights
        for black_start in list(current_black.keys()):
            for black_end in generate_moves(black_start, current_empty):
                new_black = current_black.copy()
                new_empty = current_empty.copy()
                new_black[black_end] = new_black.pop(black_start)
                new_empty.remove(black_end)
                new_empty.add(black_start)
                new_state = (frozenset(current_white.items()), frozenset(new_black.items()), frozenset(new_empty))
                if new_state not in visited:
                    visited.add(new_state)
                    queue.append((new_state, path + [f"B,{black_start},{black_end}"]))

        # Generate moves for white knights
        for white_start in list(current_white.keys()):
            for white_end in generate_moves(white_start, current_empty):
                new_white = current_white.copy()
                new_empty = current_empty.copy()
                new_white[white_end] = new_white.pop(white_start)
                new_empty.remove(white_end)
                new_empty.add(white_start)
                new_state = (frozenset(new_white.items()), frozenset(current_black.items()), frozenset(new_empty))
                if new_state not in visited:
                    visited.add(new_state)
                    queue.append((new_state, path + [f"w,{white_start},{white_end}"]))

    return "No"

# Execute the simulation
result = simulate_moves()
print(result)