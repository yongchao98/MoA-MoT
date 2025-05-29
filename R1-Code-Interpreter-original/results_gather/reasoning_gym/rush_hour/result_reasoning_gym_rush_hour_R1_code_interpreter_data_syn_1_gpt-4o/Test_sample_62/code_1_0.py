# Initial board setup
board = [
    "GBB.LM",
    "GHI.LM",
    "GHIAAN",
    "CCJK.N",
    "..JKDD",
    ".EEFF."
]

# Function to find the position of a car
def find_car(board, car):
    positions = []
    for i, row in enumerate(board):
        for j, cell in enumerate(row):
            if cell == car:
                positions.append((i, j))
    return positions

# Function to move a car
def move_car(board, car, direction, steps):
    positions = find_car(board, car)
    new_board = [list(row) for row in board]
    
    if direction == 'right':
        # Check if move is possible
        for i, j in sorted(positions, reverse=True):
            if j + steps >= len(new_board[0]) or new_board[i][j + steps] != '.':
                return board  # Return the original board if move is not possible
        # Move the car
        for i, j in sorted(positions, reverse=True):
            new_board[i][j] = '.'
            new_board[i][j + steps] = car
    elif direction == 'down':
        # Check if move is possible
        for i, j in sorted(positions, reverse=True):
            if i + steps >= len(new_board) or new_board[i + steps][j] != '.':
                return board  # Return the original board if move is not possible
        # Move the car
        for i, j in sorted(positions, reverse=True):
            new_board[i][j] = '.'
            new_board[i + steps][j] = car
    
    return [''.join(row) for row in new_board]

# Sequence of moves
moves = []

# Move N down
board = move_car(board, 'N', 'down', 1)
moves.append('N+1')

# Move M right
board = move_car(board, 'M', 'right', 1)
moves.append('M+1')

# Move L right
board = move_car(board, 'L', 'right', 1)
moves.append('L+1')

# Move AA right to exit
board = move_car(board, 'A', 'right', 2)
moves.append('A+2')

# Print the sequence of moves
print(' '.join(moves))