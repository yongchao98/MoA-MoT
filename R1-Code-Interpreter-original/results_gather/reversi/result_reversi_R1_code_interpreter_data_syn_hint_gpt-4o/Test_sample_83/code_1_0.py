def print_board(board):
    for row in board:
        print(''.join(row))

def apply_move(board, move, player):
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (-1, 1), (1, -1), (1, 1)]
    opponent = '1' if player == '0' else '0'
    x, y = move
    board[x][y] = player

    for dx, dy in directions:
        nx, ny = x + dx, y + dy
        pieces_to_flip = []

        while 0 <= nx < len(board) and 0 <= ny < len(board) and board[nx][ny] == opponent:
            pieces_to_flip.append((nx, ny))
            nx += dx
            ny += dy

        if 0 <= nx < len(board) and 0 <= ny < len(board) and board[nx][ny] == player:
            for px, py in pieces_to_flip:
                board[px][py] = player

# Initial board configuration
board = [
    ['*', '*', '*', '*', '*', '*'],
    ['*', '*', '*', '*', '0', '*'],
    ['*', '*', '1', '0', '1', '*'],
    ['*', '*', '0', '0', '*', '*'],
    ['*', '*', '*', '0', '*', '*'],
    ['*', '*', '*', '*', '*', '*']
]

# Apply the moves
apply_move(board, (4, 3), '0')  # Round 1: Your move
apply_move(board, (2, 4), '1')  # Round 2: Opponent's move
apply_move(board, (1, 4), '0')  # Round 3: Your move

# Print the final board configuration
print_board(board)