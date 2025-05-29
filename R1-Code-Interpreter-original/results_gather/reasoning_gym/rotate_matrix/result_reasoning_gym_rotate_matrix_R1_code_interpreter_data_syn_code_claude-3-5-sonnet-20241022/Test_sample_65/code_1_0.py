# Input matrix
matrix = [
    [6, 5, 0, 9, 0],
    [8, 0, 0, 4, 5],
    [0, 1, 7, 7, 1],
    [5, 3, 4, 2, 4],
    [6, 8, 3, 2, 1]
]

def rotate_270_clockwise(matrix):
    n = len(matrix)
    # For 270° clockwise, we read each column from bottom to top
    # Starting from the rightmost column
    rotated = []
    for j in range(n-1, -1, -1):  # iterate through columns from right to left
        new_row = []
        for i in range(n):  # iterate through rows from top to bottom
            new_row.append(matrix[i][j])
        rotated.append(new_row)
    return rotated

# Rotate the matrix
result = rotate_270_clockwise(matrix)

# Format output
output = "\n".join(" ".join(str(num) for num in row) for row in result)
print(f"<<<{output}>>>")