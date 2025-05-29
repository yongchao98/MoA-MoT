# Initial board setup
board = [
    ['*', '*', '*', '*'],
    ['*', '1', '0', '*'],
    ['*', '0', '0', '*'],
    ['*', '*', '0', '*']
]

# Move position
move_row, move_col = 3, 2  # (4, 3) in 1-based index

# Place the black piece
board[move_row][move_col] = '0'

# Directions: (row_delta, col_delta)
directions = [
    (-1, 0), (1, 0),  # vertical
    (0, -1), (0, 1),  # horizontal
    (-1, -1), (1, 1),  # diagonal top-left to bottom-right
    (-1, 1), (1, -1)   # diagonal top-right to bottom-left
]

# Function to check and flip pieces in a given direction
def flip_pieces(row, col, row_delta, col_delta):
    r, c = row + row_delta, col + col_delta
    pieces_to_flip = []

    while 0 <= r < 4 and 0 <= c < 4:
        if board[r][c] == '*':
            break
        if board[r][c] == '0':
            for rr, cc in pieces_to_flip:
                board[rr][cc] = '0'
            break
        pieces_to_flip.append((r, c))
        r += row_delta
        c += col_delta

# Check all directions
for dr, dc in directions:
    flip_pieces(move_row, move_col, dr, dc)

# Flatten the board for the final output
flattened_board = [cell for row in board for cell in row]
print(flattened_board)