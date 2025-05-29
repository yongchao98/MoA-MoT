# Initial matrix
matrix = [
    [0, 0, 4, 0, 0],
    [0, 0, 32, 0, 32],
    [0, 0, 0, 0, 0],
    [8, 0, 0, 0, 16],
    [0, 32, 0, 0, 4]
]

# Initial position of the number to move
position = (4, 3)

# Sequence of moves
moves = 'lluud'

# Function to perform the move
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

# Move the number in the matrix
matrix[final_position[0]][final_position[1]] = matrix[position[0]][position[1]]
matrix[position[0]][position[1]] = 0

# Print the final matrix
print(matrix)