# Initial board setup
board = [
    ".FBBB.",
    ".F.CCH",
    ".FAAGH",
    "....G.",
    "DDDEE.",
    "......"
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
    
    if direction in ['left', 'right']:
        for i, j in positions:
            new_board[i][j] = '.'
        for i, j in positions:
            if direction == 'right' and j + steps < len(new_board[0]):
                new_board[i][j + steps] = car
            elif direction == 'left' and j - steps >= 0:
                new_board[i][j - steps] = car
    elif direction in ['up', 'down']:
        for i, j in positions:
            new_board[i][j] = '.'
        for i, j in positions:
            if direction == 'down' and i + steps < len(new_board):
                new_board[i + steps][j] = car
            elif direction == 'up' and i - steps >= 0:
                new_board[i - steps][j] = car
    
    return ["".join(row) for row in new_board]

# Sequence of moves to solve the puzzle
moves = [
    ('C', 'up', 1),   # Move C up by 1
    ('H', 'right', 1), # Move H right by 1
    ('A', 'right', 2)  # Move A right by 2 to the exit
]

# Apply the moves
for car, direction, steps in moves:
    board = move_car(board, car, direction, steps)

# Print the final board state
for row in board:
    print(row)

# Print the moves in the required format
print("<<<C-1 H+1 A+2>>>")