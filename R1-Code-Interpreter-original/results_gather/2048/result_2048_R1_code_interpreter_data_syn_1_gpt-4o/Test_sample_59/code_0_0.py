# Initial matrix
matrix = [
    [32, 8, 32, 0],
    [0, 0, 32, 0],
    [0, 8, 0, 0],
    [4, 0, 0, 0]
]

# Function to move the number in the grid
def move_number(matrix, start_pos, moves):
    x, y = start_pos
    for move in moves:
        if move == 'r':
            if y < 3:
                y += 1
        elif move == 'u':
            if x > 0:
                x -= 1
        elif move == 'd':
            if x < 3:
                x += 1
        elif move == 'l':
            if y > 0:
                y -= 1
    return matrix

# Start position and moves
start_pos = (2, 2)
moves = 'rudl'

# Move the number and print the resulting matrix
result_matrix = move_number(matrix, start_pos, moves)
print(result_matrix)