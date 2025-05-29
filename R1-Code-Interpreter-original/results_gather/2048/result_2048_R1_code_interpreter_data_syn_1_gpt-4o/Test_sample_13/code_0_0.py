# Initial matrix
matrix = [
    [2, 32],
    [16, 8]
]

# Move sequence
sequence = 'ud'

# Function to move the number in the matrix
def move_number(matrix, sequence):
    # Initial position of the number to move
    pos = (0, 0)  # (row, column) for the number 2

    for move in sequence:
        if move == 'u':
            # Move up
            if pos[0] > 0:
                new_pos = (pos[0] - 1, pos[1])
                if matrix[new_pos[0]][new_pos[1]] == matrix[pos[0]][pos[1]]:
                    matrix[new_pos[0]][new_pos[1]] *= 2
                    matrix[pos[0]][pos[1]] = 0
                else:
                    matrix[new_pos[0]][new_pos[1]] = matrix[pos[0]][pos[1]]
                    matrix[pos[0]][pos[1]] = 0
                pos = new_pos
        elif move == 'd':
            # Move down
            if pos[0] < 1:
                new_pos = (pos[0] + 1, pos[1])
                if matrix[new_pos[0]][new_pos[1]] == matrix[pos[0]][pos[1]]:
                    matrix[new_pos[0]][new_pos[1]] *= 2
                    matrix[pos[0]][pos[1]] = 0
                else:
                    matrix[new_pos[0]][new_pos[1]] = matrix[pos[0]][pos[1]]
                    matrix[pos[0]][pos[1]] = 0
                pos = new_pos

    return matrix

# Apply the move sequence
result_matrix = move_number(matrix, sequence)
print(result_matrix)