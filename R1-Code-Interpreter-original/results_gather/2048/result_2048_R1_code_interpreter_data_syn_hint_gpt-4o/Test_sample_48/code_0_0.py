# Initial matrix
matrix = [
    [0, 32, 2, 2],
    [0, 0, 0, 0],
    [32, 0, 0, 0],
    [4, 4, 0, 0]
]

# Initial position of the number
pos = (2, 2)

# Move sequence
moves = 'ludl'

# Function to perform a move
def move(pos, direction):
    x, y = pos
    if direction == 'l':
        y -= 1
    elif direction == 'r':
        y += 1
    elif direction == 'u':
        x -= 1
    elif direction == 'd':
        x += 1
    return x, y

# Simulate the moves
for direction in moves:
    new_pos = move(pos, direction)
    # Check if the move is within bounds
    if 0 <= new_pos[0] < 4 and 0 <= new_pos[1] < 4:
        # Check if we can combine
        if matrix[new_pos[0]][new_pos[1]] == matrix[pos[0]][pos[1]]:
            matrix[new_pos[0]][new_pos[1]] *= 2
            matrix[pos[0]][pos[1]] = 0
        else:
            # Move the number
            matrix[new_pos[0]][new_pos[1]] = matrix[pos[0]][pos[1]]
            matrix[pos[0]][pos[1]] = 0
        pos = new_pos

# Print the final matrix
print(matrix)