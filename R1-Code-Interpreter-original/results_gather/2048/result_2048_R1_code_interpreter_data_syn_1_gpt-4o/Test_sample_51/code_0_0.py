# Initial matrix
matrix = [
    [4, 0, 0, 0],
    [0, 2, 32, 0],
    [0, 0, 0, 16],
    [0, 0, 0, 0]
]

# Initial position of the number to move
position = (3, 3)  # (4, 4) in 0-indexed

# Sequence of moves
moves = 'uuul'

# Function to move the number
def move_number(matrix, position, moves):
    x, y = position
    for move in moves:
        if move == 'u':
            x -= 1
        elif move == 'd':
            x += 1
        elif move == 'l':
            y -= 1
        elif move == 'r':
            y += 1
        # Check if the move is within bounds
        if 0 <= x < 4 and 0 <= y < 4:
            # Move the number
            matrix[x][y] = matrix[x][y] + matrix[position[0]][position[1]]
            matrix[position[0]][position[1]] = 0
            position = (x, y)
    return matrix

# Apply the moves
result_matrix = move_number(matrix, position, moves)

# Print the resulting matrix
print(result_matrix)