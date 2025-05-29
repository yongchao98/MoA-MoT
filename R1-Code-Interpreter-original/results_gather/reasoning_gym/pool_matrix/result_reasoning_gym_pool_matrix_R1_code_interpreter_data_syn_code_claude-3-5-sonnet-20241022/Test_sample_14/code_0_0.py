# Input matrix
matrix = [
    [5, 8, 4, 3, 2, 0, 2, 8],
    [0, 7, 6, 2, 2, 2, 2, 6],
    [7, 8, 8, 6, 0, 2, 4, 8]
]

# For kernel size 1, each element remains the same
# Just format to 2 decimal places
result = []
for row in matrix:
    formatted_row = [f"{x:.2f}" for x in row]
    result.append(" ".join(formatted_row))

print("\n".join(result))