# Initial matrix
matrix = [
    [4, 2, 0, 0],
    [0, 0, 0, 0],
    [2, 8, 0, 0],
    [32, 32, 4, 0]
]

# Initial position of the number
position = (3, 3)  # 0-indexed (4, 4) in 1-indexed

# Movement sequence
moves = 'dluu'

# Function to move the number
def move_number(matrix, position, moves):
    x, y = position
    for move in moves:
        if move == 'd' and x < 3:
            x += 1
        elif move == 'u' and x > 0:
            x -= 1
        elif move == 'l' and y > 0:
            y -= 1
        elif move == 'r' and y < 3:
            y += 1
    return x, y

# Get final position
final_position = move_number(matrix, position, moves)

# Update matrix
matrix[3][3] = 0  # Clear initial position
matrix[final_position[0]][final_position[1]] = 0  # Place number at final position

# Print final matrix
print(matrix)