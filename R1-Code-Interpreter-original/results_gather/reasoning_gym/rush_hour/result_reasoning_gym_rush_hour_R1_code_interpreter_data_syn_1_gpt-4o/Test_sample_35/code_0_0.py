# Initial board setup
board = [
    ['.', '.', 'B', 'B', 'B', '.'],
    ['.', '.', '.', 'F', 'G', '.'],
    ['A', 'A', '.', 'F', 'G', '.'],
    ['E', '.', '.', 'F', 'H', '.'],
    ['E', 'C', 'C', 'C', 'H', '.'],
    ['E', '.', 'D', 'D', '.', '.']
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
                    if b[i][j] == car and j + 1 < len(b[i]) and b[i][j + 1] == '.':
                        b[i][j], b[i][j + 1] = b[i][j + 1], b[i][j]
    elif direction == 'down':
        for _ in range(steps):
            for i in range(len(b) - 1, -1, -1):
                for j in range(len(b[i])):
                    if b[i][j] == car and i + 1 < len(b) and b[i + 1][j] == '.':
                        b[i][j], b[i + 1][j] = b[i + 1][j], b[i][j]

# Move F down by 1
move_car(board, 'F', 'down', 1)
# Move G down by 1
move_car(board, 'G', 'down', 1)
# Move A right by 4
move_car(board, 'A', 'right', 4)

# Print the final board
print_board(board)

# Output the moves
print("<<<F+1 G+1 A+4>>>")