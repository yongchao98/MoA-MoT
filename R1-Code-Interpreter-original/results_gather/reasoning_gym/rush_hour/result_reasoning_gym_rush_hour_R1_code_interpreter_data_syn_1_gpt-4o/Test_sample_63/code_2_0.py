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
        return False  # Car not found

    # Determine if the car is horizontal or vertical
    is_horizontal = positions[0][0] == positions[-1][0]

    # Move the car
    if is_horizontal:
        # Horizontal move
        row = positions[0][0]
        if direction == 'right':
            # Check if the path is clear
            if all(b[row][c + 1] == '.' for c in range(positions[-1][1], positions[-1][1] + steps)):
                # Move the car
                for c in range(positions[-1][1], positions[0][1] - 1, -1):
                    b[row][c + steps] = car
                    b[row][c] = '.'
                return True
        elif direction == 'left':
            # Check if the path is clear
            if all(b[row][c - 1] == '.' for c in range(positions[0][1], positions[0][1] - steps, -1)):
                # Move the car
                for c in range(positions[0][1], positions[-1][1] + 1):
                    b[row][c - steps] = car
                    b[row][c] = '.'
                return True
    else:
        # Vertical move
        col = positions[0][1]
        if direction == 'down':
            # Check if the path is clear
            if all(b[r + 1][col] == '.' for r in range(positions[-1][0], positions[-1][0] + steps)):
                # Move the car
                for r in range(positions[-1][0], positions[0][0] - 1, -1):
                    b[r + steps][col] = car
                    b[r][col] = '.'
                return True
        elif direction == 'up':
            # Check if the path is clear
            if all(b[r - 1][col] == '.' for r in range(positions[0][0], positions[0][0] - steps, -1)):
                # Move the car
                for r in range(positions[0][0], positions[-1][0] + 1):
                    b[r - steps][col] = car
                    b[r][col] = '.'
                return True

    return False

# Perform the moves
print("Initial board:")
print_board(board)

# Move M down
move_car(board, 'M', 'down', 1)
print("After moving M down:")
print_board(board)

# Move K down
move_car(board, 'K', 'down', 1)
print("After moving K down:")
print_board(board)

# Move L down
move_car(board, 'L', 'down', 1)
print("After moving L down:")
print_board(board)

# Move AA to the right
move_car(board, 'A', 'right', 3)
print("After moving AA to the right:")
print_board(board)