# Initial matrix
matrix = [
    [0, 4, 0],
    [8, 0, 2],
    [8, 2, 0]
]

# Initial position of the number
position = (1, 2)

# Move sequence
moves = 'rdl'

# Function to perform a move
def move_number(matrix, position, direction):
    x, y = position
    if direction == 'r':
        if y < 2:  # Check if we can move right
            if matrix[x][y+1] == matrix[x][y]:  # Combine if same
                matrix[x][y+1] *= 2
                matrix[x][y] = 0
            elif matrix[x][y+1] == 0:  # Move if empty
                matrix[x][y+1] = matrix[x][y]
                matrix[x][y] = 0
            y += 1
    elif direction == 'd':
        if x < 2:  # Check if we can move down
            if matrix[x+1][y] == matrix[x][y]:  # Combine if same
                matrix[x+1][y] *= 2
                matrix[x][y] = 0
            elif matrix[x+1][y] == 0:  # Move if empty
                matrix[x+1][y] = matrix[x][y]
                matrix[x][y] = 0
            x += 1
    elif direction == 'l':
        if y > 0:  # Check if we can move left
            if matrix[x][y-1] == matrix[x][y]:  # Combine if same
                matrix[x][y-1] *= 2
                matrix[x][y] = 0
            elif matrix[x][y-1] == 0:  # Move if empty
                matrix[x][y-1] = matrix[x][y]
                matrix[x][y] = 0
            y -= 1
    return (x, y)

# Apply the moves
for move in moves:
    position = move_number(matrix, position, move)

# Print the final matrix
print(matrix)