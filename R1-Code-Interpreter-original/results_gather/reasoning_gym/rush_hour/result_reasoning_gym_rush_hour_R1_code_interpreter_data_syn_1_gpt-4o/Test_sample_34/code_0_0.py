# Initial board setup
board = [
    ['G', 'B', 'B', 'B', '.', '.'],
    ['G', '.', '.', 'I', 'J', '.'],
    ['H', 'A', 'A', 'I', 'J', '.'],
    ['H', '.', '.', '.', '.', '.'],
    ['H', 'x', 'C', 'C', 'D', 'D'],
    ['E', 'E', 'F', 'F', '.', '.']
]

# Function to move a car
def move_car(board, car, direction, steps):
    # Find the car's position
    positions = [(r, c) for r in range(len(board)) for c in range(len(board[0])) if board[r][c] == car]
    if not positions:
        return False
    
    # Determine if the car is horizontal or vertical
    is_horizontal = positions[0][0] == positions[1][0]
    
    # Move the car
    if is_horizontal:
        row = positions[0][0]
        if direction == 'right':
            for _ in range(steps):
                if positions[-1][1] + 1 < len(board[0]) and board[row][positions[-1][1] + 1] == '.':
                    for r, c in reversed(positions):
                        board[r][c + 1] = car
                        board[r][c] = '.'
                    positions = [(r, c + 1) for r, c in positions]
                else:
                    return False
        elif direction == 'left':
            for _ in range(steps):
                if positions[0][1] - 1 >= 0 and board[row][positions[0][1] - 1] == '.':
                    for r, c in positions:
                        board[r][c - 1] = car
                        board[r][c] = '.'
                    positions = [(r, c - 1) for r, c in positions]
                else:
                    return False
    else:
        col = positions[0][1]
        if direction == 'down':
            for _ in range(steps):
                if positions[-1][0] + 1 < len(board) and board[positions[-1][0] + 1][col] == '.':
                    for r, c in reversed(positions):
                        board[r + 1][c] = car
                        board[r][c] = '.'
                    positions = [(r + 1, c) for r, c in positions]
                else:
                    return False
        elif direction == 'up':
            for _ in range(steps):
                if positions[0][0] - 1 >= 0 and board[positions[0][0] - 1][col] == '.':
                    for r, c in positions:
                        board[r - 1][c] = car
                        board[r][c] = '.'
                    positions = [(r - 1, c) for r, c in positions]
                else:
                    return False
    return True

# Sequence of moves
moves = [
    ('I', 'down', 1),
    ('J', 'down', 1),
    ('A', 'right', 3)
]

# Apply the moves
for car, direction, steps in moves:
    move_car(board, car, direction, steps)

# Print the final board state
for row in board:
    print(''.join(row))