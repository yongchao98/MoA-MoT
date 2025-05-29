# Initial board setup
board = [
    ['G', 'B', 'B', 'I', 'J', 'K'],
    ['G', '.', '.', 'I', 'J', 'K'],
    ['A', 'A', 'H', 'I', '.', '.'],
    ['.', '.', 'H', 'C', 'C', 'C'],
    ['.', '.', 'x', 'D', 'D', '.'],
    ['E', 'E', 'E', 'F', 'F', '.']
]

# Function to print the board
def print_board(b):
    for row in b:
        print(''.join(row))
    print()

# Function to move a car
def move_car(b, car, direction, steps):
    if direction == 'right':
        for _ in range(steps):
            for i in range(len(b)):
                for j in range(len(b[i]) - 1, -1, -1):
                    if b[i][j] == car:
                        if j + 1 < len(b[i]) and b[i][j + 1] == '.':
                            b[i][j], b[i][j + 1] = b[i][j + 1], b[i][j]
    elif direction == 'left':
        for _ in range(steps):
            for i in range(len(b)):
                for j in range(len(b[i])):
                    if b[i][j] == car:
                        if j - 1 >= 0 and b[i][j - 1] == '.':
                            b[i][j], b[i][j - 1] = b[i][j - 1], b[i][j]
    elif direction == 'down':
        for _ in range(steps):
            for i in range(len(b) - 1, -1, -1):
                for j in range(len(b[i])):
                    if b[i][j] == car:
                        if i + 1 < len(b) and b[i + 1][j] == '.':
                            b[i][j], b[i + 1][j] = b[i + 1][j], b[i][j]
    elif direction == 'up':
        for _ in range(steps):
            for i in range(len(b)):
                for j in range(len(b[i])):
                    if b[i][j] == car:
                        if i - 1 >= 0 and b[i - 1][j] == '.':
                            b[i][j], b[i - 1][j] = b[i - 1][j], b[i][j]

# Move sequence
move_car(board, 'H', 'down', 1)
move_car(board, 'C', 'right', 1)
move_car(board, 'A', 'right', 4)

# Print final board
print_board(board)