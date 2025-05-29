# Initial matrix
matrix = [
    [8, 2, 0],
    [16, 0, 0],
    [2, 0, 0]
]

# Initial position of the number to move
position = (1, 1)

# Move sequence
moves = ['r', 'u', 'r']

# Function to perform a move
def move_number(matrix, position, direction):
    x, y = position
    if direction == 'r':
        if y < 2:  # Check if not at the rightmost edge
            if matrix[x][y+1] == matrix[x][y]:  # Combine if same
                matrix[x][y+1] *= 2
                matrix[x][y] = 0
            elif matrix[x][y+1] == 0:  # Move if empty
                matrix[x][y+1] = matrix[x][y]
                matrix[x][y] = 0
            y += 1
    elif direction == 'u':
        if x > 0:  # Check if not at the top edge
            if matrix[x-1][y] == matrix[x][y]:  # Combine if same
                matrix[x-1][y] *= 2
                matrix[x][y] = 0
            elif matrix[x-1][y] == 0:  # Move if empty
                matrix[x-1][y] = matrix[x][y]
                matrix[x][y] = 0
            x -= 1
    return (x, y)

# Execute the moves
for move in moves:
    position = move_number(matrix, position, move)

# Print the final matrix
print(matrix)