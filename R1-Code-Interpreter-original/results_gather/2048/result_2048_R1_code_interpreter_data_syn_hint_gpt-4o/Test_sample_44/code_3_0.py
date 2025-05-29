# Initial matrix
matrix = [
    [2, 4, 32, 0],
    [0, 0, 16, 0],
    [0, 0, 0, 0],
    [2, 8, 0, 0]
]

# Move sequence
moves = 'lrud'

def move_left(matrix):
    for i in range(4):
        new_row = [num for num in matrix[i] if num != 0]
        for j in range(1, len(new_row)):
            if new_row[j] == new_row[j - 1]:
                new_row[j - 1] *= 2
                new_row[j] = 0
        new_row = [num for num in new_row if num != 0]
        matrix[i] = new_row + [0] * (4 - len(new_row))

def move_right(matrix):
    for i in range(4):
        new_row = [num for num in matrix[i] if num != 0]
        for j in range(len(new_row) - 2, -1, -1):
            if new_row[j] == new_row[j + 1]:
                new_row[j + 1] *= 2
                new_row[j] = 0
        new_row = [num for num in new_row if num != 0]
        matrix[i] = [0] * (4 - len(new_row)) + new_row

def move_up(matrix):
    for j in range(4):
        new_col = [matrix[i][j] for i in range(4) if matrix[i][j] != 0]
        for i in range(1, len(new_col)):
            if new_col[i] == new_col[i - 1]:
                new_col[i - 1] *= 2
                new_col[i] = 0
        new_col = [num for num in new_col if num != 0]
        for i in range(4):
            matrix[i][j] = new_col[i] if i < len(new_col) else 0

def move_down(matrix):
    for j in range(4):
        new_col = [matrix[i][j] for i in range(4) if matrix[i][j] != 0]
        for i in range(len(new_col) - 2, -1, -1):
            if new_col[i] == new_col[i + 1]:
                new_col[i + 1] *= 2
                new_col[i] = 0
        new_col = [num for num in new_col if num != 0]
        for i in range(4):
            matrix[i][j] = new_col[i - (4 - len(new_col))] if i >= 4 - len(new_col) else 0

# Simulate the moves
for move in moves:
    if move == 'l':
        move_left(matrix)
    elif move == 'r':
        move_right(matrix)
    elif move == 'u':
        move_up(matrix)
    elif move == 'd':
        move_down(matrix)

# Print the final matrix
print(matrix)