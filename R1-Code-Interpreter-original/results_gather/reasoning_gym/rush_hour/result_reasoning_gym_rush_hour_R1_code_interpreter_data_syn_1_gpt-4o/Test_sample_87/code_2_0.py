# Initial board setup
board = [
    ['G', '.', '.', 'I', 'B', 'B'],
    ['G', '.', 'H', 'I', '.', 'J'],
    ['G', '.', 'H', 'A', 'A', 'J'],
    ['C', 'C', 'D', 'D', '.', 'K'],
    ['E', 'E', 'E', 'F', 'F', 'K'],
    ['.', '.', '.', '.', '.', '.']
]

# Function to move a car
def move_car(board, car, direction, steps):
    # Find the car's position
    positions = [(i, j) for i in range(6) for j in range(6) if board[i][j] == car]
    if not positions:
        return False

    # Determine if the car is horizontal or vertical
    is_horizontal = positions[0][0] == positions[1][0]

    # Calculate new positions
    if is_horizontal:
        if direction == 'right':
            new_positions = [(i, j + steps) for i, j in positions]
        else:
            new_positions = [(i, j - steps) for i, j in positions]
    else:
        if direction == 'down':
            new_positions = [(i + steps, j) for i, j in positions]
        else:
            new_positions = [(i - steps, j) for i, j in positions]

    # Check if new positions are within bounds and empty
    for i, j in new_positions:
        if i < 0 or i >= 6 or j < 0 or j >= 6 or (board[i][j] != '.' and board[i][j] != car):
            return False

    # Move the car
    for i, j in positions:
        board[i][j] = '.'

    for i, j in new_positions:
        board[i][j] = car

    return True

# Sequence of moves
moves = []

# Move J down
if move_car(board, 'J', 'down', 1):
    moves.append('J+1')

# Move AA to the right
if move_car(board, 'A', 'right', 2):
    moves.append('A+2')

# Print the sequence of moves
print(' '.join(moves))