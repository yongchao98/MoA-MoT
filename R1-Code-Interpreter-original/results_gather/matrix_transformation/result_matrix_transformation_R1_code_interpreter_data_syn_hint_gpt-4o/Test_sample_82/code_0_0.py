def rotate_90_clockwise(matrix):
    n = len(matrix)
    rotated_matrix = []
    for j in range(n):
        new_row = []
        for i in range(n-1, -1, -1):
            new_row.append(matrix[i][j])
        rotated_matrix.append(new_row)
    return rotated_matrix

matrix = [
    ['*', '&', '&', '@'],
    ['&', '%', '*', '!'],
    ['!', '*', '%', '^'],
    ['^', '@', '%', '$']
]

rotated_matrix = rotate_90_clockwise(matrix)
print(rotated_matrix)