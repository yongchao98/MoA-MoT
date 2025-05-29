# Initial matrix
matrix = [
    [4, 8],
    [4, 2]
]

# Initial position of the number to move
position = (1, 0)  # (row, column) in zero-indexed

# Sequence of moves
moves = 'dd'

# Function to move the number in the matrix
def move_number(matrix, position, moves):
    for move in moves:
        if move == 'd':
            # Move down
            if position[0] < len(matrix) - 1:
                # Swap the numbers if they are the same
                if matrix[position[0]][position[1]] == matrix[position[0] + 1][position[1]]:
                    matrix[position[0] + 1][position[1]] *= 2
                    matrix[position[0]][position[1]] = 0
                else:
                    # Just move the number down
                    matrix[position[0] + 1][position[1]], matrix[position[0]][position[1]] = matrix[position[0]][position[1]], matrix[position[0] + 1][position[1]]
                position = (position[0] + 1, position[1])
    return matrix

# Apply the moves
result_matrix = move_number(matrix, position, moves)

# Print the resulting matrix
print(result_matrix)