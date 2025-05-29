# Initial matrix
matrix = [
    [0, 16, 0, 32, 0],
    [0, 0, 0, 0, 0],
    [32, 2, 8, 32, 0],
    [0, 32, 0, 0, 0],
    [0, 0, 0, 0, 0]
]

# Initial position of the number
position = (3, 4)  # 0-indexed position (4, 5)

# Sequence of moves
moves = 'rdrrl'

# Function to move the number in the matrix
def move_number(matrix, position, moves):
    x, y = position
    for move in moves:
        if move == 'r' and y < 4:
            y += 1
        elif move == 'd' and x < 4:
            x += 1
        elif move == 'l' and y > 0:
            y -= 1
        elif move == 'u' and x > 0:
            x -= 1
    return x, y

# Get the final position after all moves
final_position = move_number(matrix, position, moves)

# Since the number is 0, the matrix remains unchanged
# Print the final matrix
print(matrix)