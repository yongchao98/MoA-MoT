# Input matrix
matrix = [
    [0, 5, 6],
    [9, 4, 0],
    [0, 4, 0],
    [4, 1, 6],
    [9, 9, 3],
    [9, 4, 6],
    [6, 9, 4]
]

# For kernel size 1, we just format the output
rows = len(matrix)
cols = len(matrix[0])

# Format and print the result
result = ""
for i in range(rows):
    row = " ".join(f"{x:.2f}" for x in matrix[i])
    result += row + "\n"

print(result.strip())