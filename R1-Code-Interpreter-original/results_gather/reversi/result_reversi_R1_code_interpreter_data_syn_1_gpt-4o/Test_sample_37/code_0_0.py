# Initial board setup
board = [
    ['*', '*', '*', '*'],
    ['*', '1', '0', '*'],
    ['*', '1', '0', '0'],
    ['*', '1', '*', '*']
]

# Function to flip pieces in a given direction
def flip_pieces(board, start, direction, player):
    x, y = start
    dx, dy = direction
    to_flip = []
    while 0 <= x < 4 and 0 <= y < 4:
        x += dx
        y += dy
        if not (0 <= x < 4 and 0 <= y < 4):
            return
        if board[x][y] == '*':
            return
        if board[x][y] == player:
            for fx, fy in to_flip:
                board[fx][fy] = player
            return
        to_flip.append((x, y))

# Directions: right, left, down, up, diagonals
directions = [(0, 1), (0, -1), (1, 0), (-1, 0), (1, 1), (1, -1), (-1, 1), (-1, -1)]

# Round 1: Your move (black) at (3,4)
board[2][3] = '0'
for direction in directions:
    flip_pieces(board, (2, 3), direction, '0')

# Round 2: Opponent's move (white) at (4,2)
board[3][1] = '1'
for direction in directions:
    flip_pieces(board, (3, 1), direction, '1')

# Flatten the board to the required format
flattened_board = [cell for row in board for cell in row]
print(flattened_board)