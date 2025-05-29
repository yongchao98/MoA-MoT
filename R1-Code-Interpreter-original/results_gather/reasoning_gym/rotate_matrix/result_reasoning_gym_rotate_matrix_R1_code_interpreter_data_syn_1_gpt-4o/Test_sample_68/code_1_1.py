# Define the matrix
matrix = [
    [7, 7, 2, 6],
    [5, 6, 8, 1],
    [8, 1, 1, 7],
    [4, 0, 3, 5]
]

# Transpose the matrix
transposed_matrix = list(zip(*matrix))

# Reverse the order of the rows to achieve 90-degree counterclockwise rotation
rotated_matrix = [list(row)[::-1] for row in transposed_matrix]

# Print the rotated matrix
print(rotated_matrix)