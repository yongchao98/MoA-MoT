# Input matrix
matrix = [
    [4, 0, 1, 5, 0, 3],
    [1, 2, 7, 0, 3, 2]
]

kernel_size = 3
rows = len(matrix)
cols = len(matrix[0])

# For kernel size 3, we'll have 2 regions to average
result = []

# First region (first 3 columns)
elements = [
    matrix[0][0], matrix[0][1], matrix[0][2],
    matrix[1][0], matrix[1][1], matrix[1][2]
]
avg1 = sum(elements) / len(elements)

# Second region (last 3 columns)
elements = [
    matrix[0][3], matrix[0][4], matrix[0][5],
    matrix[1][3], matrix[1][4], matrix[1][5]
]
avg2 = sum(elements) / len(elements)

print(f"{avg1:.2f} {avg2:.2f}")