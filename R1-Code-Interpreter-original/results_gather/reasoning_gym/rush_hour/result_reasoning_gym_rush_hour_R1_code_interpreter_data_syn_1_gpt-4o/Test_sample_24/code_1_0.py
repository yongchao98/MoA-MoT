# Initial board setup
board = [
    ['.', 'E', 'B', 'B', 'B', '.'],
    ['.', 'E', '.', '.', 'F', '.'],
    ['A', 'A', '.', '.', 'F', '.'],
    ['.', '.', 'C', 'C', 'C', '.'],
    ['.', '.', 'D', 'D', 'D', '.'],
    ['.', '.', '.', '.', '.', '.']
]

# Function to print the board
def print_board(board):
    for row in board:
        print(''.join(row))
    print()

# Function to move a car
def move_car(board, car, direction, steps):
    # Find the car's position
    positions = [(r, c) for r in range(6) for c in range(6) if board[r][c] == car]
    if not positions:
        return False

    # Determine the orientation of the car
    vertical = all(p[1] == positions[0][1] for p in positions)
    horizontal = all(p[0] == positions[0][0] for p in positions)

    # Calculate new positions
    new_positions = []
    if vertical:
        for r, c in positions:
            new_r = r + steps if direction == 'down' else r - steps
            new_positions.append((new_r, c))
    elif horizontal:
        for r, c in positions:
            new_c = c + steps if direction == 'right' else c - steps
            new_positions.append((r, new_c))

    # Check if new positions are valid
    for r, c in new_positions:
        if r < 0 or r >= 6 or c < 0 or c >= 6 or (board[r][c] != '.' and board[r][c] != car):
            return False

    # Move the car
    for r, c in positions:
        board[r][c] = '.'
    for r, c in new_positions:
        board[r][c] = car

    return True

# Sequence of moves
moves = [
    ('F', 'down', 1),  # Move F down by 1
    ('C', 'up', 1),    # Move C up by 1
    ('A', 'right', 4)  # Move A right by 4 to the exit
]

# Apply the moves
for car, direction, steps in moves:
    if not move_car(board, car, direction, steps):
        print(f"Failed to move car {car} {direction} by {steps} steps.")
    print_board(board)

# Final board state
print("Final board state:")
print_board(board)