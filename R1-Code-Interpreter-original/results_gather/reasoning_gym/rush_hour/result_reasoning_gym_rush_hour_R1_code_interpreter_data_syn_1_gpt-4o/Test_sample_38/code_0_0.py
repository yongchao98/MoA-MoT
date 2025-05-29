# Initial board setup
board = [
    ['G', '.', 'J', 'B', 'B', '.'],
    ['G', '.', 'J', 'K', '.', 'L'],
    ['H', 'A', 'A', 'K', '.', 'L'],
    ['H', 'I', 'C', 'C', '.', 'x'],
    ['H', 'I', 'D', 'D', '.', '.'],
    ['E', 'E', 'E', 'F', 'F', '.']
]

# Function to print the board
def print_board(b):
    for row in b:
        print(''.join(row))
    print()

# Function to move a car
def move_car(b, car, direction, spaces):
    # Find the car's position
    positions = [(r, c) for r in range(len(b)) for c in range(len(b[0])) if b[r][c] == car]
    if not positions:
        return False  # Car not found

    # Determine if the car is horizontal or vertical
    is_horizontal = positions[0][0] == positions[1][0]

    # Move the car
    if is_horizontal:
        # Horizontal car
        if direction == '+':
            # Move right
            for r, c in reversed(positions):
                if c + spaces >= len(b[0]) or b[r][c + spaces] not in ('.', car):
                    return False  # Invalid move
                b[r][c + spaces] = car
                b[r][c] = '.'
        elif direction == '-':
            # Move left
            for r, c in positions:
                if c - spaces < 0 or b[r][c - spaces] not in ('.', car):
                    return False  # Invalid move
                b[r][c - spaces] = car
                b[r][c] = '.'
    else:
        # Vertical car
        if direction == '+':
            # Move down
            for r, c in reversed(positions):
                if r + spaces >= len(b) or b[r + spaces][c] not in ('.', car):
                    return False  # Invalid move
                b[r + spaces][c] = car
                b[r][c] = '.'
        elif direction == '-':
            # Move up
            for r, c in positions:
                if r - spaces < 0 or b[r - spaces][c] not in ('.', car):
                    return False  # Invalid move
                b[r - spaces][c] = car
                b[r][c] = '.'

    return True

# Simulate the moves
moves = [('K', '+', 1), ('L', '+', 1), ('A', '+', 3)]
for car, direction, spaces in moves:
    if not move_car(board, car, direction, spaces):
        print(f"Move {car}{direction}{spaces} is invalid.")
    print_board(board)

# Output the moves
print("<<<K+1 L+1 A+3>>>")