# Initial matrix
matrix = [
    [0, 2, 0, 0, 4],
    [32, 16, 0, 2, 32],
    [0, 2, 0, 0, 0],
    [0, 2, 0, 0, 0],
    [8, 0, 2, 0, 0]
]

# Starting position
x, y = 2, 1

# Move sequence
moves = 'ddldlu'

# Directions mapping
directions = {
    'u': (-1, 0),
    'd': (1, 0),
    'l': (0, -1),
    'r': (0, 1)
}

# Function to perform the move
def move_number(matrix, x, y, moves):
    for move in moves:
        dx, dy = directions[move]
        nx, ny = x + dx, y + dy
        
        # Check if the new position is within bounds
        if 0 <= nx < 5 and 0 <= ny < 5:
            if matrix[nx][ny] == matrix[x][y]:
                # Combine numbers
                matrix[nx][ny] *= 2
                matrix[x][y] = 0
            elif matrix[nx][ny] == 0:
                # Move number
                matrix[nx][ny] = matrix[x][y]
                matrix[x][y] = 0
            x, y = nx, ny
    return matrix

# Perform the moves
result_matrix = move_number(matrix, x, y, moves)

# Print the resulting matrix
print(result_matrix)