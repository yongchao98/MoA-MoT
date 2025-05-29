# Initial matrix
matrix = [
    [0, 0, 0, 0, 2],
    [0, 0, 0, 0, 0],
    [0, 0, 8, 0, 32],
    [0, 0, 0, 0, 4],
    [4, 32, 0, 0, 0]
]

# Starting position
x, y = 2, 2

# Sequence of moves
moves = 'rdlluu'

# Move directions
directions = {
    'r': (0, 1),
    'd': (1, 0),
    'l': (0, -1),
    'u': (-1, 0)
}

# Simulate the moves
for move in moves:
    dx, dy = directions[move]
    new_x, new_y = x + dx, y + dy
    
    # Check if the new position is within bounds
    if 0 <= new_x < 5 and 0 <= new_y < 5:
        # Check if the numbers can combine
        if matrix[new_x][new_y] == matrix[x][y]:
            matrix[new_x][new_y] *= 2
            matrix[x][y] = 0
        elif matrix[new_x][new_y] == 0:
            matrix[new_x][new_y] = matrix[x][y]
            matrix[x][y] = 0
        
        # Update the position
        x, y = new_x, new_y

# Print the final matrix
print(matrix)