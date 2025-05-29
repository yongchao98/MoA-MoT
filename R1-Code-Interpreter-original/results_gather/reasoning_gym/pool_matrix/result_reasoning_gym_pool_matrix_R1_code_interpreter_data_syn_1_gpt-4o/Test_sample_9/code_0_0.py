# Define the input matrix
input_matrix = [
    [0, 3, 2, 0, 1, 3, 9],
    [7, 5, 8, 3, 3, 2, 6],
    [0, 4, 6, 9, 3, 8, 4],
    [7, 6, 1, 5, 2, 1, 5],
    [3, 2, 4, 1, 2, 3, 4],
    [5, 4, 5, 7, 6, 2, 4],
    [5, 8, 4, 6, 1, 4, 9],
    [9, 5, 2, 4, 4, 4, 3]
]

# Kernel size
kernel_size = 1

# Perform average pooling
output_matrix = []
for i in range(0, len(input_matrix), kernel_size):
    row = []
    for j in range(0, len(input_matrix[0]), kernel_size):
        # Since kernel size is 1, each element is its own average
        row.append(round(input_matrix[i][j], 2))
    output_matrix.append(row)

# Print the output matrix
print(output_matrix)