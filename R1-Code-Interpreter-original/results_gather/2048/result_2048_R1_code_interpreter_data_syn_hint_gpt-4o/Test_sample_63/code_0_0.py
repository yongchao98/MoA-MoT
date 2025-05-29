# Initial matrix
matrix = [
    [2, 0, 0, 0, 0],
    [0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0],
    [2, 2, 2, 8, 0]
]

# Initial position of the number 8
position = (4, 3)

# Sequence of moves
moves = 'ldrll'

# Function to move the number in the matrix
def move_number(matrix, position, moves):
    x, y = position
    for move in moves:
        if move == 'l':
            if y > 0:
                y -= 1
        elif move == 'r':
            if y < 4:
                y += 1
        elif move == 'u':
            if x > 0:
                x -= 1
        elif move == 'd':
            if x < 4:
                x += 1
    return x, y

# Get the final position after all moves
final_position = move_number(matrix, position, moves)

# Update the matrix with the new position of the number 8
matrix[4][3] = 0  # Clear the initial position
matrix[final_position[0]][final_position[1]] = 8  # Set the new position

# Print the final matrix
print(matrix)