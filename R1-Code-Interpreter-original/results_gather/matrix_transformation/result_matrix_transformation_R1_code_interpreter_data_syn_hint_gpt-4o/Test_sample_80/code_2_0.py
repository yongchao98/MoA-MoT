# Original matrix
matrix = [
    ['^', '&', '@', '%'],
    ['*', '!', '!', '#'],
    ['!', '#', '*', '%'],
    ['*', '$', '^', '@']
]

# Number of rows (or columns, since it's a square matrix)
n = len(matrix)

# Initialize the rotated matrix
rotated_matrix = [['' for _ in range(n)] for _ in range(n)]

# Perform the rotation
for i in range(n):
    for j in range(n):
        rotated_matrix[j][n-1-i] = matrix[i][j]

# Flatten the rotated matrix for the output format
flattened_result = [item for sublist in rotated_matrix for item in sublist]

# Format the result as a string with the required format
formatted_result = "<<<" + ",".join(flattened_result) + ">>>"

# Print the formatted result
print(formatted_result)