# Initial matrix
matrix = [
    [0, 0, 16, 0],
    [2, 8, 0, 0],
    [0, 4, 4, 0],
    [0, 0, 0, 0]
]

# Initial position of the number to move
position = (3, 0)  # (4, 1) in 0-indexed

# Sequence of moves
moves = 'uuul'

# Function to move the number in the matrix
def move_number(matrix, position, moves):
    x, y = position
    for move in moves:
        if move == 'u' and x > 0:
            x -= 1
        elif move == 'd' and x < len(matrix) - 1:
            x += 1
        elif move == 'l' and y > 0:
            y -= 1
        elif move == 'r' and y < len(matrix[0]) - 1:
            y += 1
    return x, y

# Get the final position after all moves
final_position = move_number(matrix, position, moves)

# Move the number in the matrix
matrix[final_position[0]][final_position[1]] = matrix[position[0]][position[1]]
matrix[position[0]][position[1]] = 0

# Print the final matrix
print(matrix)