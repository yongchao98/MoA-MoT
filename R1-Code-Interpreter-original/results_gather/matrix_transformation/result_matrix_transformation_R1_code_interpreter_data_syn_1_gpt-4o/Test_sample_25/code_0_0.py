# Original matrix
matrix = [
    ["6", "7", "5"],
    ["9", "2", "4"]
]

# Number of rows and columns
rows = len(matrix)
cols = len(matrix[0])

# Create a new matrix for the rotated result
rotated_matrix = [["" for _ in range(rows)] for _ in range(cols)]

# Perform the rotation
for i in range(rows):
    for j in range(cols):
        rotated_matrix[j][rows - 1 - i] = matrix[i][j]

# Flatten the rotated matrix and format the output
flattened_result = [item for sublist in rotated_matrix for item in sublist]
print(flattened_result)