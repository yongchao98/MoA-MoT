from collections import deque

# Define the initial board state
initial_board = [
    ['.', '.', 'w', '.'],
    ['.', 'w', '.', '.'],
    ['.', 'B', 'B', '.']
]

# Define the target board state
target_board = [
    ['.', 'B', 'B', '.'],
    ['.', '.', '.', '.'],
    ['.', '.', 'w', 'w']
]

# Define knight moves (L-shape)
knight_moves = [
    (2, 1), (2, -1), (-2, 1), (-2, -1),
    (1, 2), (1, -2), (-1, 2), (-1, -2)
]

# Convert board position to a string for easy comparison
def board_to_string(board):
    return ''.join(''.join(row) for row in board)

# Find the positions of all knights
def find_knights(board):
    white_knights = []
    black_knights = []
    for r in range(3):
        for c in range(4):
            if board[r][c] == 'w':
                white_knights.append((r, c))
            elif board[r][c] == 'B':
                black_knights.append((r, c))
    return white_knights, black_knights

# Check if a position is within the board
def is_within_board(r, c):
    return 0 <= r < 3 and 0 <= c < 4

# Perform BFS to find the sequence of moves
def bfs_solve(initial_board, target_board):
    initial_state = board_to_string(initial_board)
    target_state = board_to_string(target_board)
    queue = deque([(initial_board, [], 'w')])  # (current board, moves, next turn)
    visited = set()
    visited.add(initial_state)

    while queue:
        current_board, moves, turn = queue.popleft()
        current_state = board_to_string(current_board)

        if current_state == target_state:
            return moves

        # Find current knights
        white_knights, black_knights = find_knights(current_board)
        knights = white_knights if turn == 'w' else black_knights

        # Try all possible moves for the current turn
        for r, c in knights:
            for dr, dc in knight_moves:
                nr, nc = r + dr, c + dc
                if is_within_board(nr, nc) and current_board[nr][nc] == '.':
                    # Make the move
                    new_board = [row[:] for row in current_board]
                    new_board[nr][nc] = turn
                    new_board[r][c] = '.'
                    new_state = board_to_string(new_board)
                    if new_state not in visited:
                        visited.add(new_state)
                        new_moves = moves + [f"{turn},{chr(c + 65)}{3 - r},{chr(nc + 65)}{3 - nr}"]
                        queue.append((new_board, new_moves, 'B' if turn == 'w' else 'w'))

    return "No"

# Solve the puzzle
result = bfs_solve(initial_board, target_board)
print(result)