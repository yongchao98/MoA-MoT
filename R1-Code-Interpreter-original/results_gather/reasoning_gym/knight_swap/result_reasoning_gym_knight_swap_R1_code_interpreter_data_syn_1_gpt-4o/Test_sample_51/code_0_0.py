from collections import deque

# Define the board and initial positions
board = {
    'A2': 'w', 'C2': 'w',
    'D3': 'B', 'C1': 'B'
}

# Define target positions
target_positions = {
    'w': ['D3', 'C1'],
    'B': ['A2', 'C2']
}

# Define knight moves
knight_moves = [
    (2, 1), (2, -1), (-2, 1), (-2, -1),
    (1, 2), (1, -2), (-1, 2), (-1, -2)
]

# Convert board positions to coordinates
def pos_to_coord(pos):
    col, row = pos
    return ord(col) - ord('A'), int(row) - 1

# Convert coordinates to board positions
def coord_to_pos(coord):
    col, row = coord
    return chr(col + ord('A')) + str(row + 1)

# Check if a position is within the board
def is_within_board(coord):
    col, row = coord
    return 0 <= col < 4 and 0 <= row < 3

# BFS to find the sequence of moves
def bfs_knight_swap():
    queue = deque([(board, [], 'w')])  # (current board, moves, current turn)
    visited = set()

    while queue:
        current_board, moves, turn = queue.popleft()
        board_tuple = tuple(sorted(current_board.items()))
        
        if board_tuple in visited:
            continue
        visited.add(board_tuple)

        # Check if the current board matches the target positions
        if all(current_board.get(pos) == 'w' for pos in target_positions['w']) and \
           all(current_board.get(pos) == 'B' for pos in target_positions['B']):
            return moves

        # Generate possible moves for the current turn
        for pos, knight in current_board.items():
            if knight != turn:
                continue

            current_coord = pos_to_coord(pos)
            for move in knight_moves:
                new_coord = (current_coord[0] + move[0], current_coord[1] + move[1])
                new_pos = coord_to_pos(new_coord)

                if is_within_board(new_coord) and new_pos not in current_board:
                    # Make the move
                    new_board = current_board.copy()
                    del new_board[pos]
                    new_board[new_pos] = knight

                    # Add the move to the queue
                    new_moves = moves + [f"{knight},{pos},{new_pos}"]
                    next_turn = 'B' if turn == 'w' else 'w'
                    queue.append((new_board, new_moves, next_turn))

    return "No"

# Find the solution
solution = bfs_knight_swap()
print(solution)