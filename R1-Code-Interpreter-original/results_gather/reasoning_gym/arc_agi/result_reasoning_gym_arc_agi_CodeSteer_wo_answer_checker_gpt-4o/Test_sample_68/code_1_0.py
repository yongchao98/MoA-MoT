import numpy as np

def rotate_90_clockwise(matrix):
    np_matrix = np.array(matrix)
    # Transpose the matrix and then reverse each row to achieve a 90-degree clockwise rotation
    rotated_matrix = np.rot90(np_matrix, k=-1)
    return rotated_matrix

input_grid = [
    [1, 3, 5],
    [3, 5, 3],
    [8, 8, 3]
]

output_grid = rotate_90_clockwise(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))