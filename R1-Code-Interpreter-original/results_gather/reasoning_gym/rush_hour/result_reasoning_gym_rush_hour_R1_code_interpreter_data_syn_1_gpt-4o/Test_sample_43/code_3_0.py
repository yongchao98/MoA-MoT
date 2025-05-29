# Initial board setup
board = [
    ['B', 'B', 'B', '.', 'x', 'K'],
    ['C', 'C', 'D', 'D', 'D', 'K'],
    ['H', 'A', 'A', 'J', '.', '.'],
    ['H', '.', '.', 'J', 'E', 'E'],
    ['.', '.', 'I', 'F', 'F', 'x'],
    ['G', 'G', 'I', '.', '.', '.']
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

    # Determine the orientation of the car
    if len(set(r for r, c in positions)) == 1:
        orientation = 'H'  # Horizontal
    else:
        orientation = 'V'  # Vertical

    # Move the car
    if orientation == 'H':
        if direction == '+':
            # Move right
            for r, c in sorted(positions, reverse=True):
                if c + steps < 6 and b[r][c + steps] == '.':
                    b[r][c + steps] = car
                    b[r][c] = '.'
                else:
                    return False
        elif direction == '-':
            # Move left
            for r, c in sorted(positions):
                if c - steps >= 0 and b[r][c - steps] == '.':
                    b[r][c - steps] = car
                    b[r][c] = '.'
                else:
                    return False
    elif orientation == 'V':
        if direction == '+':
            # Move down
            for r, c in sorted(positions, reverse=True):
                if r + steps < 6 and b[r + steps][c] == '.':
                    b[r + steps][c] = car
                    b[r][c] = '.'
                else:
                    return False
        elif direction == '-':
            # Move up
            for r, c in sorted(positions):
                if r - steps >= 0 and b[r - steps][c] == '.':
                    b[r - steps][c] = car
                    b[r][c] = '.'
                else:
                    return False
    return True

# Sequence of moves
moves = [
    ('K', '+', 1),  # Move K down by 1 space
    ('J', '+', 1),  # Move J down by 1 space
    ('D', '+', 1),  # Move D right
    ('A', '+', 2)   # Move A right to exit
]

# Apply the moves
for car, direction, steps in moves:
    if not move_car(board, car, direction, steps):
        print("Invalid move for car", car)
        break
    print_board(board)

# Final sequence of moves
print("<<<K+1 J+1 D+1 A+2>>>")