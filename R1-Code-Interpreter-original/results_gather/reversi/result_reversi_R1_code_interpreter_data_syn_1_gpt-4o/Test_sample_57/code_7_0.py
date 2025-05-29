# Initial board configuration
board = [
    ['*', '0', '*', '*'],
    ['*', '0', '0', '*'],
    ['1', '0', '1', '*'],
    ['0', '*', '*', '*']
]

# Function to flip pieces in a line
def flip_line(board, start, end, player):
    x1, y1 = start
    x2, y2 = end
    dx = (x2 - x1) // max(abs(x2 - x1), 1)
    dy = (y2 - y1) // max(abs(y2 - y1), 1)
    x, y = x1 + dx, y1 + dy
    while (x, y) != (x2, y2):
        board[x][y] = player
        x += dx
        y += dy

# Function to apply a move
def apply_move(board, move, player):
    x, y = move
    board[x][y] = player
    opponent = '1' if player == '0' else '0'
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (1, 1), (-1, 1), (1, -1)]
    for dx, dy in directions:
        nx, ny = x + dx, y + dy
        if 0 <= nx < 4 and 0 <= ny < 4 and board[nx][ny] == opponent:
            # Check if there's a player's piece in the line
            line = []
            while 0 <= nx < 4 and 0 <= ny < 4:
                if board[nx][ny] == player:
                    flip_line(board, (x, y), (nx, ny), player)
                    break
                elif board[nx][ny] == '*':
                    break
                nx += dx
                ny += dy

# Apply the moves
apply_move(board, (0, 1), '0')  # Round 1
apply_move(board, (2, 0), '1')  # Round 2
apply_move(board, (3, 0), '0')  # Round 3

# Flatten the board for the final output
flattened_board = [cell for row in board for cell in row]
print(flattened_board)