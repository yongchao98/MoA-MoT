# Initial matrix
matrix = [
    [0, 4, 8, 0],
    [0, 2, 0, 0],
    [0, 32, 2, 0],
    [2, 4, 0, 0]
]

# Initial position of the number to move
position = (2, 2)

# Sequence of moves
moves = 'rrud'

# Function to perform the moves
def move_number(matrix, position, moves):
    x, y = position
    for move in moves:
        if move == 'r':
            y += 1
        elif move == 'l':
            y -= 1
        elif move == 'u':
            x -= 1
        elif move == 'd':
            x += 1
        
        # Check if the current position is within bounds and if the cell is empty
        if 0 <= x < 4 and 0 <= y < 4 and matrix[x][y] == 0:
            # Move the number
            matrix[x][y] = matrix[position[0]][position[1]]
            # Clear the previous position
            matrix[position[0]][position[1]] = 0
            # Update the current position
            position = (x, y)
    
    return matrix

# Perform the moves
result_matrix = move_number(matrix, position, moves)

# Print the resulting matrix
print(result_matrix)