# Define the matrix
matrix = [
    [5, 0, 1, 8, 4, 6, 1, 8, 1, 6],
    [7, 0, 9, 3, 3, 5, 0, 2, 6, 4]
]

# Identity transformation (no change)
identity_matrix = matrix

# Mirror along the counterdiagonal
# For a 2x10 matrix, mirroring along the counterdiagonal swaps elements across the line from (0,9) to (1,0)
mirrored_matrix = [
    [matrix[1][9], matrix[1][8], matrix[1][7], matrix[1][6], matrix[1][5], matrix[1][4], matrix[1][3], matrix[1][2], matrix[1][1], matrix[1][0]],
    [matrix[0][9], matrix[0][8], matrix[0][7], matrix[0][6], matrix[0][5], matrix[0][4], matrix[0][3], matrix[0][2], matrix[0][1], matrix[0][0]]
]

# Rotate 360 degrees (no change)
final_matrix = mirrored_matrix

# Print the final matrix
for row in final_matrix:
    print(' '.join(map(str, row)))