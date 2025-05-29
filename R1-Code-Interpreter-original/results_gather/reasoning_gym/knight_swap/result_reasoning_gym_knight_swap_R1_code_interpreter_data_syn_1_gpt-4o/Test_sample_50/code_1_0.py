from collections import deque

# Define the initial positions of the knights
initial_white_knights = ('A2', 'C1')
initial_black_knights = ('D3', 'D2')
initial_empty_squares = ('B3', 'B2')

# Define the possible moves for a knight
knight_moves = [(2, 1), (2, -1), (-2, 1), (-2, -1), (1, 2), (1, -2), (-1, 2), (-1, -2)]

# Convert board positions to coordinates
def pos_to_coord(pos):
    col, row = pos
    return ord(col) - ord('A'), int(row) - 1

# Convert coordinates back to board positions
def coord_to_pos(coord):
    col, row = coord
    return chr(col + ord('A')) + str(row + 1)

# Check if a move is valid
def is_valid_move(start, end, occupied):
    if end in occupied:
        return False
    start_coord = pos_to_coord(start)
    end_coord = pos_to_coord(end)
    for move in knight_moves:
        if (start_coord[0] + move[0], start_coord[1] + move[1]) == end_coord:
            return True
    return False

# Perform BFS to find the sequence of moves
def bfs_knight_swap():
    # Initial state
    initial_state = (initial_white_knights, initial_black_knights, initial_empty_squares)
    queue = deque([(initial_state, [])])
    visited = set()
    visited.add(initial_state)

    while queue:
        (white_knights, black_knights, empty_squares), moves = queue.popleft()

        # Check if the swap is complete
        if white_knights == initial_black_knights and black_knights == initial_white_knights:
            return moves

        # Generate all possible moves for white knights
        for i, white_pos in enumerate(white_knights):
            for empty_pos in empty_squares:
                if is_valid_move(white_pos, empty_pos, white_knights + black_knights):
                    new_white_knights = list(white_knights)
                    new_empty_squares = list(empty_squares)
                    new_white_knights[i] = empty_pos
                    new_empty_squares[new_empty_squares.index(empty_pos)] = white_pos
                    new_state = (tuple(new_white_knights), black_knights, tuple(new_empty_squares))
                    if new_state not in visited:
                        visited.add(new_state)
                        queue.append((new_state, moves + [f'w,{white_pos},{empty_pos}']))

        # Generate all possible moves for black knights
        for i, black_pos in enumerate(black_knights):
            for empty_pos in empty_squares:
                if is_valid_move(black_pos, empty_pos, white_knights + black_knights):
                    new_black_knights = list(black_knights)
                    new_empty_squares = list(empty_squares)
                    new_black_knights[i] = empty_pos
                    new_empty_squares[new_empty_squares.index(empty_pos)] = black_pos
                    new_state = (white_knights, tuple(new_black_knights), tuple(new_empty_squares))
                    if new_state not in visited:
                        visited.add(new_state)
                        queue.append((new_state, moves + [f'B,{black_pos},{empty_pos}']))

    return "No"

# Execute the BFS search
result = bfs_knight_swap()
print(result)