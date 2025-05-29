# Initial board state
board = [
    ['*', '*', '*', '*'],
    ['*', '1', '0', '*'],
    ['*', '0', '0', '*'],
    ['*', '*', '0', '*']
]

# Place the black piece at (4, 3)
board[3][2] = '0'

# Directions: (row_change, col_change)
directions = [
    (0, 1), (0, -1),  # Horizontal
    (1, 0), (-1, 0),  # Vertical
    (1, 1), (-1, -1), # Diagonal \
    (1, -1), (-1, 1)  # Diagonal /
]

# Function to flip pieces in a given direction
def flip_pieces(r, c, dr, dc):
    path = []
    r += dr
    c += dc
    while 0 <= r < 4 and 0 <= c < 4:
        if board[r][c] == '*':
            return
        if board[r][c] == '0':
            for pr, pc in path:
                board[pr][pc] = '0'
            return
        path.append((r, c))
        r += dr
        c += dc

# Check all directions from the new piece
for dr, dc in directions:
    flip_pieces(3, 2, dr, dc)

# Flatten the board to a single list
flattened_board = [cell for row in board for cell in row]
print(flattened_board)