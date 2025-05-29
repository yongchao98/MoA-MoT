# Initial matrix
matrix = [
    [2, 4, 32, 0],
    [0, 0, 16, 0],
    [0, 0, 0, 0],
    [2, 8, 0, 0]
]

# Move sequence
moves = 'lrud'

def compress_and_merge(line):
    """Compress non-zero elements and merge equal adjacent elements."""
    new_line = [num for num in line if num != 0]
    for i in range(len(new_line) - 1):
        if new_line[i] == new_line[i + 1]:
            new_line[i] *= 2
            new_line[i + 1] = 0
    return [num for num in new_line if num != 0]

def move_left(matrix):
    for i in range(4):
        new_row = compress_and_merge(matrix[i])
        matrix[i] = new_row + [0] * (4 - len(new_row))

def move_right(matrix):
    for i in range(4):
        new_row = compress_and_merge(matrix[i][::-1])
        matrix[i] = [0] * (4 - len(new_row)) + new_row[::-1]

def move_up(matrix):
    for j in range(4):
        new_col = compress_and_merge([matrix[i][j] for i in range(4)])
        for i in range(4):
            matrix[i][j] = new_col[i] if i < len(new_col) else 0

def move_down(matrix):
    for j in range(4):
        new_col = compress_and_merge([matrix[i][j] for i in range(4)][::-1])
        for i in range(4):
            matrix[i][j] = new_col[::-1][i] if i < len(new_col) else 0

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