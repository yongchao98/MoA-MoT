# Original matrix
matrix = [
    ["7", "1", "2"],
    ["7", "8", "9"]
]

# Number of rows and columns
rows = len(matrix)
cols = len(matrix[0])

# Initialize the rotated matrix with empty strings
rotated_matrix = [["" for _ in range(rows)] for _ in range(cols)]

# Fill the rotated matrix
for i in range(rows):
    for j in range(cols):
        rotated_matrix[j][rows - 1 - i] = matrix[i][j]

# Flatten the rotated matrix to a list of strings
flattened_result = [item for sublist in rotated_matrix for item in sublist]

# Print the result in the required format
print(f"<<<{','.join(flattened_result)}>>>")