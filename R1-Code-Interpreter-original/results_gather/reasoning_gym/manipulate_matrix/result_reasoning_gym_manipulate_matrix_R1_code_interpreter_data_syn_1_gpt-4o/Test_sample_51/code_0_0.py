# Initial matrix
matrix = [
    [1, 2, 6, 0, 0, 3, 3, 2, 4],
    [1, 0, 6, 2, 5, 4, 8, 3, 7],
    [3, 7, 2, 0, 0, 0, 7, 5, 6],
    [9, 2, 8, 5, 7, 7, 9, 5, 7],
    [3, 2, 6, 8, 8, 2, 9, 2, 0],
    [5, 7, 3, 7, 9, 0, 1, 7, 9],
    [8, 6, 2, 3, 4, 4, 1, 5, 1],
    [2, 1, 9, 6, 9, 9, 1, 2, 4],
    [7, 4, 5, 6, 7, 2, 9, 5, 9]
]

# Step 2: Map each occurrence of 4 to 0
for i in range(len(matrix)):
    for j in range(len(matrix[i])):
        if matrix[i][j] == 4:
            matrix[i][j] = 0

# Step 3: Map each occurrence of 7 to 6
for i in range(len(matrix)):
    for j in range(len(matrix[i])):
        if matrix[i][j] == 7:
            matrix[i][j] = 6

# Step 4: Crop the matrix to rows 2-2 and columns 8-8 (1-indexed)
# This corresponds to row index 1 and column index 7 in 0-indexed Python
cropped_matrix = [[matrix[1][7]]]

# Step 5: Map each occurrence of 7 to 2 in the cropped matrix
for i in range(len(cropped_matrix)):
    for j in range(len(cropped_matrix[i])):
        if cropped_matrix[i][j] == 7:
            cropped_matrix[i][j] = 2

# Output the final matrix
print(cropped_matrix)