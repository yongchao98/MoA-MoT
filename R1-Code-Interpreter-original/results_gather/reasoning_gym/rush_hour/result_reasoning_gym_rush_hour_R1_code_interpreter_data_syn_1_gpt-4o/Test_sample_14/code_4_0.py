# Initial board setup
board = [
    ['.', '.', 'I', 'B', 'B', '.'],
    ['.', 'G', 'I', '.', 'K', '.'],
    ['.', 'G', 'A', 'A', 'K', '.'],
    ['F', 'G', 'J', 'C', 'C', '.'],
    ['F', 'H', 'J', '.', '.', '.'],
    ['F', 'H', 'D', 'D', 'E', 'E']
]

# Function to print the board
def print_board(b):
    for row in b:
        print(''.join(row))
    print()

# Function to move a car
def move_car(b, car, direction, steps):
    # Find the car's position
    positions = [(r, c) for r in range(6) for c in range(6) if b[r][c] == car]
    if not positions:
        return False

    # Determine orientation
    if len(set(r for r, c in positions)) == 1:
        orientation = 'H'  # Horizontal
    else:
        orientation = 'V'  # Vertical

    # Calculate new positions
    if orientation == 'H':
        if direction == '+':
            new_positions = [(r, c + steps) for r, c in positions]
        else:
            new_positions = [(r, c - steps) for r, c in positions]
    else:
        if direction == '+':
            new_positions = [(r + steps, c) for r, c in positions]
        else:
            new_positions = [(r - steps, c) for r, c in positions]

    # Check if new positions are valid
    for r, c in new_positions:
        if r < 0 or r >= 6 or c < 0 or c >= 6 or (b[r][c] != '.' and b[r][c] != car):
            print(f"Invalid move for car {car} to position ({r}, {c})")
            return False

    # Move the car
    for r, c in positions:
        b[r][c] = '.'
    for r, c in new_positions:
        b[r][c] = car

    return True

# Sequence of moves to solve the puzzle
moves = [
    ('K', '+', 1),  # Move K right
    ('C', '-', 1),  # Move C left
    ('J', '-', 1),  # Move J left
    ('A', '+', 2)   # Move A right to the exit
]

# Apply the moves
for car, direction, steps in moves:
    success = move_car(board, car, direction, steps)
    print(f"Move {car} {direction}{steps}: {'Success' if success else 'Failed'}")
    print_board(board)

# Print the final board
print("Final board:")
print_board(board)