# Define the initial positions of the knights
white_knights = {'B3': 'w', 'A1': 'w'}
black_knights = {'D2': 'B', 'C1': 'B'}
empty_squares = {'A2', 'D3', 'B1'}

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

# Check if a move is valid
def is_valid_move(start, end, empty_squares):
    start_coord = pos_to_coord(start)
    end_coord = pos_to_coord(end)
    if end in empty_squares:
        for move in knight_moves:
            if (start_coord[0] + move[0], start_coord[1] + move[1]) == end_coord:
                return True
    return False

# Check if a position is within the board
def is_within_board(coord):
    col, row = coord
    return 0 <= col < 4 and 0 <= row < 3

# Simulate the moves
def simulate_moves():
    moves = []
    current_empty = empty_squares.copy()
    current_white = white_knights.copy()
    current_black = black_knights.copy()

    # Black moves first
    for black_start in list(current_black.keys()):
        for move in knight_moves:
            black_end_coord = (pos_to_coord(black_start)[0] + move[0], pos_to_coord(black_start)[1] + move[1])
            if is_within_board(black_end_coord):
                black_end = coord_to_pos(black_end_coord)
                if is_valid_move(black_start, black_end, current_empty):
                    moves.append(f"B,{black_start},{black_end}")
                    current_empty.remove(black_end)
                    current_empty.add(black_start)
                    current_black[black_end] = current_black.pop(black_start)
                    break

    # Alternate moves between white and black
    while current_white != target_black or current_black != target_white:
        for white_start in list(current_white.keys()):
            for move in knight_moves:
                white_end_coord = (pos_to_coord(white_start)[0] + move[0], pos_to_coord(white_start)[1] + move[1])
                if is_within_board(white_end_coord):
                    white_end = coord_to_pos(white_end_coord)
                    if is_valid_move(white_start, white_end, current_empty):
                        moves.append(f"w,{white_start},{white_end}")
                        current_empty.remove(white_end)
                        current_empty.add(white_start)
                        current_white[white_end] = current_white.pop(white_start)
                        break

        for black_start in list(current_black.keys()):
            for move in knight_moves:
                black_end_coord = (pos_to_coord(black_start)[0] + move[0], pos_to_coord(black_start)[1] + move[1])
                if is_within_board(black_end_coord):
                    black_end = coord_to_pos(black_end_coord)
                    if is_valid_move(black_start, black_end, current_empty):
                        moves.append(f"B,{black_start},{black_end}")
                        current_empty.remove(black_end)
                        current_empty.add(black_start)
                        current_black[black_end] = current_black.pop(black_start)
                        break

        if current_white == target_black and current_black == target_white:
            return moves

    return "No"

# Execute the simulation
result = simulate_moves()
print(result)