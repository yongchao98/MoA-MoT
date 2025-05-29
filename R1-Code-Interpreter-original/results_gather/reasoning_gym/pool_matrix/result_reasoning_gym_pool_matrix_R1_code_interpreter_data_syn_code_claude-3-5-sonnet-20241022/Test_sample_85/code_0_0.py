# Input matrix
matrix = [
    [0, 0, 1, 7],
    [3, 1, 5, 1],
    [0, 2, 9, 2],
    [6, 2, 5, 1],
    [1, 6, 7, 6]
]

# For kernel size 1, each element remains the same
# Just format to 2 decimal places
result = []
for row in matrix:
    formatted_row = [f"{x:.2f}" for x in row]
    result.append(" ".join(formatted_row))

print("\n".join(result))