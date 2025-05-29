def print_board(board):
    for row in board:
        print(''.join(row))

def apply_move(board, move, player):
    n = len(board)
    x, y = move
    board[x][y] = player

    # Directions: horizontal, vertical, and diagonal
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (1, 1), (-1, 1), (1, -1)]

    for dx, dy in directions:
        nx, ny = x + dx, y + dy
        pieces_to_flip = []

        while 0 <= nx < n and 0 <= ny < n:
            if board[nx][ny] == '*':
                break
            if board[nx][ny] == player:
                for px, py in pieces_to_flip:
                    board[px][py] = player
                break
            pieces_to_flip.append((nx, ny))
            nx += dx
            ny += dy

# Initial board configuration
board = [
    ['*', '*', '*', '0'],
    ['*', '1', '0', '1'],
    ['*', '0', '0', '*'],
    ['*', '*', '0', '*']
]

# Moves: (row, column) are 0-indexed
moves = [(3, 2, '0'), (1, 3, '1'), (0, 3, '0')]

# Apply each move
for move in moves:
    x, y, player = move
    apply_move(board, (x, y), player)

# Print the final board configuration
print_board(board)