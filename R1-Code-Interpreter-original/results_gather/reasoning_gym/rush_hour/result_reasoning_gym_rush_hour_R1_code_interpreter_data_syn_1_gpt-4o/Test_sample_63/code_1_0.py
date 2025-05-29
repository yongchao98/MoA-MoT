# Initial board setup
board = [
    ['.', 'B', 'B', 'B', 'K', 'M'],
    ['C', 'C', '.', 'I', 'K', 'M'],
    ['A', 'A', '.', 'I', 'L', 'M'],
    ['G', 'D', 'D', 'J', 'L', '.'],
    ['G', '.', 'H', 'J', 'E', 'E'],
    ['F', 'F', 'H', '.', '.', '.']
]

# Function to print the board
def print_board(b):
    for row in b:
        print(''.join(row))
    print()

# Function to move a car
def move_car(b, car, direction, steps):
    # Find the car's position
    positions = [(r, c) for r in range(len(b)) for c in range(len(b[0])) if b[r][c] == car]
    if not positions:
        return False

    # Determine if the car is horizontal or vertical
    is_horizontal = positions[0][0] == positions[1][0]

    # Move the car
    if is_horizontal:
        if direction == 'right':
            for r, c in reversed(positions):
                if c + steps >= len(b[0]) or b[r][c + steps] not in ('.', car):
                    return False
                b[r][c], b[r][c + steps] = '.', car
        elif direction == 'left':
            for r, c in positions:
                if c - steps < 0 or b[r][c - steps] not in ('.', car):
                    return False
                b[r][c], b[r][c - steps] = '.', car
    else:
        if direction == 'down':
            for r, c in reversed(positions):
                if r + steps >= len(b) or b[r + steps][c] not in ('.', car):
                    return False
                b[r][c], b[r + steps][c] = '.', car
        elif direction == 'up':
            for r, c in positions:
                if r - steps < 0 or b[r - steps][c] not in ('.', car):
                    return False
                b[r][c], b[r - steps][c] = '.', car

    return True

# Simulate the moves
print("Initial board:")
print_board(board)

# Move M down by 1
move_car(board, 'M', 'down', 1)
print("After moving M down by 1:")
print_board(board)

# Move AA right by 3
move_car(board, 'A', 'right', 3)
print("After moving AA right by 3:")
print_board(board)

# Output the moves
print("<<<M+1 A+3>>>")