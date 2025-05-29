# Initial board setup
board = [
    ['.', '.', 'B', 'B', 'B', 'I'],
    ['.', '.', '.', '.', 'H', 'I'],
    ['A', 'A', '.', '.', 'H', 'I'],
    ['C', 'C', 'C', '.', 'H', '.'],
    ['.', 'G', 'D', 'D', 'E', 'E'],
    ['.', 'G', '.', 'F', 'F', '.']
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
    if is_horizontal:
        for r, c in positions:
            new_c = c + spaces if direction == 'right' else c - spaces
            if new_c < 0 or new_c >= len(b[0]) or (b[r][new_c] != '.' and b[r][new_c] != car):
                return False  # Invalid move
            new_positions.append((r, new_c))
    else:
        for r, c in positions:
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

# Simulate the moves
print("Initial board:")
print_board(board)

# Move H down by 1
move_car(board, 'H', 'down', 1)
print("After moving H down by 1:")
print_board(board)

# Move I down by 1
move_car(board, 'I', 'down', 1)
print("After moving I down by 1:")
print_board(board)

# Move A right by 3
move_car(board, 'A', 'right', 3)
print("After moving A right by 3:")
print_board(board)