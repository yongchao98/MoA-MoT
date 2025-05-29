# Input matrix
matrix = [
    [0, 3, 2, 0, 1, 3, 9],
    [7, 5, 8, 3, 3, 2, 6],
    [0, 4, 6, 9, 3, 8, 4],
    [7, 6, 1, 5, 2, 1, 5],
    [3, 2, 4, 1, 2, 3, 4],
    [5, 4, 5, 7, 6, 2, 4],
    [5, 8, 4, 6, 1, 4, 9],
    [9, 5, 2, 4, 4, 4, 3]
]

# For kernel size 1, average pooling will keep the same values
# but we need to format them to 2 decimal places
result = []
for row in matrix:
    formatted_row = [f"{x:.2f}" for x in row]
    print(" ".join(formatted_row))