# Initial matrix
matrix = [
    [2, 0, 8, 2],
    [0, 16, 0, 0],
    [0, 0, 8, 4],
    [0, 16, 32, 0]
]

# Function to move right
def move_right(row):
    new_row = [x for x in row if x != 0]  # Remove zeros
    for i in range(len(new_row) - 1, 0, -1):
        if new_row[i] == new_row[i - 1]:
            new_row[i] *= 2
            new_row[i - 1] = 0
    new_row = [x for x in new_row if x != 0]  # Remove zeros again
    return [0] * (len(row) - len(new_row)) + new_row

# Function to move left
def move_left(row):
    new_row = [x for x in row if x != 0]  # Remove zeros
    for i in range(len(new_row) - 1):
        if new_row[i] == new_row[i + 1]:
            new_row[i] *= 2
            new_row[i + 1] = 0
    new_row = [x for x in new_row if x != 0]  # Remove zeros again
    return new_row + [0] * (len(row) - len(new_row))

# Apply the sequence 'rrll' to the third row
matrix[2] = move_right(matrix[2])
matrix[2] = move_right(matrix[2])
matrix[2] = move_left(matrix[2])
matrix[2] = move_left(matrix[2])

# Print the final matrix
print(matrix)