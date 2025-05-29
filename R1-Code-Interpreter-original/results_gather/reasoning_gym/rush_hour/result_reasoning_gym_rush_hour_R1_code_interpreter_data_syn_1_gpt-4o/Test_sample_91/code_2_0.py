# Initial board setup
board = [
    ['B', 'B', 'C', 'C', 'C', '.'],
    ['.', '.', '.', '.', 'H', 'I'],
    ['.', '.', 'A', 'A', 'H', 'I'],
    ['.', '.', '.', 'G', 'D', 'D'],
    ['.', '.', '.', 'G', '.', '.'],
    ['x', 'E', 'E', 'F', 'F', '.']
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

    # Calculate new positions
    new_positions = []
    for r, c in positions:
        if is_horizontal:
            new_c = c + spaces if direction == 'right' else c - spaces
            if new_c < 0 or new_c >= len(b[0]) or (b[r][new_c] != '.' and b[r][new_c] != car):
                return False  # Invalid move
            new_positions.append((r, new_c))
        else:
            new_r = r + spaces if direction == 'down' else r - spaces
            if new_r < 0 or new_r >= len(b) or (b[new_r][c] != '.' and b[new_r][c] != car):
                return False  # Invalid move
            new_positions.append((new_r, c))

    # Move the car
    for r, c in positions:
        b[r][c] = '.'
    for r, c in new_positions:
        b[r][c] = car

    return True

# Perform the moves
moves = [
    ('H', 'up', 1),  # Move 'HI' up
    ('G', 'up', 1),  # Move 'G' up
    ('D', 'right', 1), # Move 'DD' right
    ('A', 'right', 2)  # Move 'AA' to the exit
]

# Apply the moves
for car, direction, spaces in moves:
    if not move_car(board, car, direction, spaces):
        print(f"Failed to move {car} {direction} by {spaces} spaces.")
    print_board(board)

# Final board state
print_board(board)