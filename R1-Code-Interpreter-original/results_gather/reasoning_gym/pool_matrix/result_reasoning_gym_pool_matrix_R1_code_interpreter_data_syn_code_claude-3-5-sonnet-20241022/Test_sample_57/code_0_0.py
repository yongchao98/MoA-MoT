# Input matrix
matrix = [
    [3, 4, 6, 5, 8, 2],
    [1, 7, 8, 0, 7, 6],
    [4, 1, 6, 8, 6, 6]
]

# For kernel size 1, we just need to format each number to 2 decimal places
result = []
for row in matrix:
    formatted_row = [f"{x:.2f}" for x in row]
    result.append(" ".join(formatted_row))

print("\n".join(result))