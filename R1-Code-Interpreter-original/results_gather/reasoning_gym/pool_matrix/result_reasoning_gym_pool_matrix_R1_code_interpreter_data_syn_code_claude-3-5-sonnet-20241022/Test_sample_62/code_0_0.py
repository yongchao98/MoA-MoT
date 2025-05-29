# Input matrix
matrix = [
    [2, 0, 5, 3],
    [2, 6, 9, 8],
    [3, 6, 0, 8],
    [7, 0, 2, 7],
    [9, 1, 8, 7],
    [0, 8, 8, 8],
    [3, 6, 2, 7],
    [2, 5, 0, 1],
    [0, 1, 7, 3],
    [1, 2, 4, 2]
]

# For kernel size 1, output is same as input but formatted to 2 decimal places
output = []
for row in matrix:
    formatted_row = [f"{x:.2f}" for x in row]
    print(" ".join(formatted_row))