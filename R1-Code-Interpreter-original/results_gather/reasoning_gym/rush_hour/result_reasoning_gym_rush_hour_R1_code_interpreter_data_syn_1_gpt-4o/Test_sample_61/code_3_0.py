# Initial board setup
board = [
    ['B', 'B', 'B', 'C', 'C', 'M'],
    ['D', 'D', 'E', 'E', 'E', 'M'],
    ['A', 'A', 'J', 'K', '.', '.'],
    ['I', '.', 'J', 'K', 'F', 'F'],
    ['I', '.', 'G', 'G', 'L', 'x'],
    ['H', 'H', '.', '.', 'L', '.']
]

# Function to print the board
def print_board(b):
    for row in b:
        print(''.join(row))
    print()

# Function to move a car
def move_car(b, car, direction, spaces):
    if direction == 'right':
        for _ in range(spaces):
            for i in range(len(b)):
                for j in range(len(b[i]) - 1, -1, -1):
                    if b[i][j] == car and j + 1 < len(b[i]) and b[i][j + 1] == '.':
                        b[i][j], b[i][j + 1] = b[i][j + 1], b[i][j]
    elif direction == 'down':
        for _ in range(spaces):
            for i in range(len(b) - 1, -1, -1):
                for j in range(len(b[i])):
                    if b[i][j] == car and i + 1 < len(b) and b[i + 1][j] == '.':
                        b[i][j], b[i + 1][j] = b[i + 1][j], b[i][j]

# Simulate the moves
print("Initial Board:")
print_board(board)

# Move M down by 1
move_car(board, 'M', 'down', 1)
# Move E down by 1
move_car(board, 'E', 'down', 1)
# Move K down by 1
move_car(board, 'K', 'down', 1)
# Move AA to the right by 4
move_car(board, 'A', 'right', 4)

print("Final Board:")
print_board(board)