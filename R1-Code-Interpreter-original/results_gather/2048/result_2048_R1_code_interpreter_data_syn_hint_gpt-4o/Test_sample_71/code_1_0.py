# Initial matrix
matrix = [
    [0, 0, 0, 0, 0],
    [2, 0, 0, 0, 0],
    [0, 0, 0, 0, 32],
    [0, 4, 0, 0, 16],
    [16, 0, 0, 32, 0]
]

# Initial position of the number to move
x, y = 4, 3

# Sequence of moves
moves = 'ldrll'

# Function to perform a move
def move(x, y, direction):
    if direction == 'l':
        return x, max(0, y - 1)
    elif direction == 'r':
        return x, min(4, y + 1)
    elif direction == 'u':
        return max(0, x - 1), y
    elif direction == 'd':
        return min(4, x + 1), y
    return x, y

# Simulate the moves
for move_direction in moves:
    new_x, new_y = move(x, y, move_direction)
    if matrix[new_x][new_y] == matrix[x][y]:
        # Combine numbers if they are the same
        matrix[new_x][new_y] *= 2
        matrix[x][y] = 0
    elif matrix[new_x][new_y] == 0:
        # Move the number if the new position is empty
        matrix[new_x][new_y] = matrix[x][y]
        matrix[x][y] = 0
    x, y = new_x, new_y

# Print the final matrix
print(matrix)