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
    x, y = move(x, y, move_direction)

# Update the matrix with the final position
matrix[4][3] = 0  # Clear the initial position
matrix[x][y] = 0  # Place the number at the new position

# Print the final matrix
print(matrix)