# Input matrix
matrix = [
    [9, 2, 5],
    [8, 9, 5],
    [4, 9, 2]
]

def rotate_90_clockwise(mat):
    N = len(mat)
    # Create new matrix with reversed rows
    return [[mat[N-1-j][i] for j in range(N)] for i in range(N)]

# Since 450 degrees = 90 degrees (as 450 = 360 + 90)
rotated = rotate_90_clockwise(matrix)

# Format output
output = '\n'.join(' '.join(str(x) for x in row) for row in rotated)
print(f"<<<{output}>>>")