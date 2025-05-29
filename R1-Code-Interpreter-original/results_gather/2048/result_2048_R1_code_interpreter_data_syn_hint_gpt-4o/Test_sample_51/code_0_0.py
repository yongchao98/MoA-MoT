# Initial matrix
matrix = [
    [4, 0, 0, 0],
    [0, 2, 32, 0],
    [0, 0, 0, 16],
    [0, 0, 0, 0]
]

# Sequence of moves: 'uuul'
moves = ['u', 'u', 'u', 'l']

# Initial position of the number to move
position = (3, 3)  # (4, 4) in 0-indexed

# Function to move the number
def move_number(matrix, position, moves):
    for move in moves:
        x, y = position
        if move == 'u' and x > 0:
            position = (x - 1, y)
        elif move == 'd' and x < 3:
            position = (x + 1, y)
        elif move == 'l' and y > 0:
            position = (x, y - 1)
        elif move == 'r' and y < 3:
            position = (x, y + 1)
    return matrix

# Apply the moves
result_matrix = move_number(matrix, position, moves)

# Print the resulting matrix
print(result_matrix)