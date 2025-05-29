# Initial board setup
board = [
    ['E', '.', '.', 'H', '.', '.'],
    ['E', '.', '.', 'H', 'x', 'x'],
    ['E', 'A', 'A', 'I', 'J', '.'],
    ['.', 'F', 'G', 'I', 'J', 'K'],
    ['.', 'F', 'G', 'B', 'B', 'K'],
    ['C', 'C', 'D', 'D', 'D', '.']
]

# Function to move a car
def move_car(board, car, direction, spaces):
    # Find the car's position
    positions = [(r, c) for r in range(len(board)) for c in range(len(board[0])) if board[r][c] == car]
    if not positions:
        return False  # Car not found

    # Determine if the car is horizontal or vertical
    is_horizontal = positions[0][0] == positions[1][0]

    # Move the car
    if is_horizontal:
        if direction == 'right':
            for r, c in reversed(positions):
                if c + spaces >= len(board[0]) or board[r][c + spaces] not in ('.', car):
                    return False  # Move is invalid
                board[r][c + spaces] = car
                board[r][c] = '.'
        elif direction == 'left':
            for r, c in positions:
                if c - spaces < 0 or board[r][c - spaces] not in ('.', car):
                    return False  # Move is invalid
                board[r][c - spaces] = car
                board[r][c] = '.'
    else:
        if direction == 'down':
            for r, c in reversed(positions):
                if r + spaces >= len(board) or board[r + spaces][c] not in ('.', car):
                    return False  # Move is invalid
                board[r + spaces][c] = car
                board[r][c] = '.'
        elif direction == 'up':
            for r, c in positions:
                if r - spaces < 0 or board[r - spaces][c] not in ('.', car):
                    return False  # Move is invalid
                board[r - spaces][c] = car
                board[r][c] = '.'

    return True

# Perform the moves
moves = [
    ('I', 'up', 1),
    ('J', 'left', 1),
    ('K', 'down', 1),
    ('A', 'right', 3)
]

for car, direction, spaces in moves:
    if not move_car(board, car, direction, spaces):
        print("Invalid move for car", car)
        break

# Print the final board state
for row in board:
    print(''.join(row))