# Initial matrix
matrix = [
    [8, 0, 4],
    [32, 16, 16],
    [0, 0, 2]
]

# Move sequence
moves = 'lll'

# Function to move left
def move_left(matrix, row):
    # Create a new row to store the result of the move
    new_row = [0, 0, 0]
    index = 0
    last_value = 0
    
    for col in range(3):
        if matrix[row][col] != 0:
            if last_value == 0:
                last_value = matrix[row][col]
            elif last_value == matrix[row][col]:
                new_row[index] = last_value * 2
                index += 1
                last_value = 0
            else:
                new_row[index] = last_value
                index += 1
                last_value = matrix[row][col]
    
    if last_value != 0:
        new_row[index] = last_value
    
    # Update the matrix row
    matrix[row] = new_row

# Apply the moves
for move in moves:
    if move == 'l':
        move_left(matrix, 2)

# Print the final matrix
print(matrix)