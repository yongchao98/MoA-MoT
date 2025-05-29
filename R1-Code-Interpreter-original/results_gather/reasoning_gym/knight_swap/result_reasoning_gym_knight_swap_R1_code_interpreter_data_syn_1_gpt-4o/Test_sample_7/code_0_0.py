from collections import deque

# Define the board and initial positions
board = {
    'A1': None, 'B1': None, 'C1': 'B', 'D1': None,
    'A2': None, 'B2': None, 'C2': 'w', 'D2': 'w',
    'A3': None, 'B3': 'B', 'C3': None, 'D3': None
}

# Define knight moves
knight_moves = [
    (2, 1), (2, -1), (-2, 1), (-2, -1),
    (1, 2), (1, -2), (-1, 2), (-1, -2)
]

# Convert board position to coordinates
def pos_to_coord(pos):
    col, row = pos
    return ord(col) - ord('A'), int(row) - 1

# Convert coordinates to board position
def coord_to_pos(x, y):
    return chr(x + ord('A')) + str(y + 1)

# Check if a move is valid
def is_valid_move(x, y):
    return 0 <= x < 4 and 0 <= y < 3

# BFS to find the sequence of moves
def bfs_solve():
    queue = deque()
    queue.append((board, [], 'w'))  # (current board, moves, current turn)
    visited = set()

    while queue:
        current_board, moves, turn = queue.popleft()
        board_tuple = tuple(current_board.items())
        
        if board_tuple in visited:
            continue
        visited.add(board_tuple)

        # Check if the current board is the target configuration
        if current_board['B3'] == 'w' and current_board['C1'] == 'w' and current_board['C2'] == 'B' and current_board['D2'] == 'B':
            return moves

        # Find all knights of the current turn
        for pos, knight in current_board.items():
            if knight == turn:
                x, y = pos_to_coord(pos)
                # Try all possible knight moves
                for dx, dy in knight_moves:
                    nx, ny = x + dx, y + dy
                    if is_valid_move(nx, ny):
                        new_pos = coord_to_pos(nx, ny)
                        if current_board[new_pos] is None:  # Move to empty square
                            new_board = current_board.copy()
                            new_board[new_pos] = turn
                            new_board[pos] = None
                            new_moves = moves + [f"{turn},{pos},{new_pos}"]
                            next_turn = 'B' if turn == 'w' else 'w'
                            queue.append((new_board, new_moves, next_turn))
    
    return "No"

# Find the solution
solution = bfs_solve()
print(solution)