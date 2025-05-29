# Original matrix
matrix = [
    ["6", "P", "T", "D"],
    ["7", "L", "L", "4"],
    ["P", "C", "O", "S"]
]

# Number of rows and columns
n = len(matrix)
m = len(matrix[0])

# Initialize the rotated matrix
rotated_matrix = []

# Perform the rotation
for j in range(m):
    new_row = []
    for i in range(n-1, -1, -1):
        new_row.append(matrix[i][j])
    rotated_matrix.append(new_row)

# Flatten the rotated matrix for output
flattened_result = [item for sublist in rotated_matrix for item in sublist]

# Print the result
print(flattened_result)