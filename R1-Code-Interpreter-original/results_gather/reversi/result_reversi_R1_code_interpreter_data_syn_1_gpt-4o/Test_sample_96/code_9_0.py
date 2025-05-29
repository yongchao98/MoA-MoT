# Initial board configuration
board = [
    ['*', '*', '*', '*', '*', '*'],
    ['*', '*', '0', '*', '*', '*'],
    ['*', '*', '0', '1', '1', '*'],
    ['*', '*', '0', '0', '*', '*'],
    ['*', '*', '*', '0', '*', '*'],
    ['*', '*', '*', '*', '*', '*']
]

# Function to flip pieces in a given direction
def flip_pieces(board, x, y, dx, dy, player):
    opponent = '1' if player == '0' else '0'
    to_flip = []
    i, j = x + dx, y + dy
    while 0 <= i < 6 and 0 <= j < 6:
        if board[i][j] == opponent:
            to_flip.append((i, j))
        elif board[i][j] == player:
            for fx, fy in to_flip:
                board[fx][fy] = player
            break
        else:
            break
        i += dx
        j += dy

# Function to apply a move
def apply_move(board, x, y, player):
    board[x][y] = player
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (1, 1), (-1, 1), (1, -1)]
    for dx, dy in directions:
        flip_pieces(board, x, y, dx, dy, player)

# Apply the moves
apply_move(board, 4, 3, '0')  # Round 1: Your move
apply_move(board, 2, 4, '1')  # Round 2: Opponent's move
apply_move(board, 1, 2, '0')  # Round 3: Your move

# Print the final board configuration
final_board = [cell for row in board for cell in row]
print(final_board)