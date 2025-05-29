# Input matrix
matrix = [
    [6, 8, 6, 5, 5, 8],
    [6, 5, 0, 8, 6, 3],
    [9, 0, 6, 3, 2, 8],
    [9, 0, 2, 7, 9, 2],
    [2, 3, 5, 5, 9, 2],
    [3, 4, 4, 1, 5, 1]
]

# For 900 degrees (which is equivalent to 180 degrees),
# we just need to reverse the matrix both row-wise and column-wise
rotated = matrix[::-1]  # reverse rows
rotated = [row[::-1] for row in rotated]  # reverse each row

# Format the output
output = ""
for row in rotated:
    output += " ".join(map(str, row)) + "\n"
print(output.strip())