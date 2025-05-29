# Initial matrix
matrix = [
    [8, 2, 0],
    [16, 0, 0],
    [2, 0, 0]
]

# Move sequence
moves = 'rur'

# Function to move the number
def move_number(matrix, start_pos, moves):
    x, y = start_pos
    for move in moves:
        if move == 'r' and y < 2:
            y += 1
        elif move == 'l' and y > 0:
            y -= 1
        elif move == 'u' and x > 0:
            x -= 1
        elif move == 'd' and x < 2:
            x += 1
    return x, y

# Initial position of the number to move
start_pos = (1, 0)

# Get the final position after moves
final_pos = move_number(matrix, start_pos, moves)

# Move the number in the matrix
number_to_move = matrix[start_pos[0]][start_pos[1]]
matrix[start_pos[0]][start_pos[1]] = 0
matrix[final_pos[0]][final_pos[1]] = number_to_move

# Print the final matrix
print(matrix)