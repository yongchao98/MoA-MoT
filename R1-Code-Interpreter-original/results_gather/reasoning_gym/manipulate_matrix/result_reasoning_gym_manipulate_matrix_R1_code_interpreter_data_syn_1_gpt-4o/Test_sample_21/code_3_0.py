# Initial matrix
matrix = [
    [8, 7, 8, 5, 7, 0, 0, 2, 5],
    [2, 4, 0, 2, 0, 5, 3, 7, 3],
    [3, 2, 2, 5, 8, 5, 8, 6, 0],
    [5, 4, 2, 1, 8, 6, 6, 7, 3],
    [2, 0, 7, 7, 3, 5, 9, 3, 2]
]

# Rotate the matrix 90 degrees clockwise
rotated_matrix = list(zip(*matrix[::-1]))

# Remove every 2-nd row (1-indexed)
final_matrix = [row for i, row in enumerate(rotated_matrix) if (i + 1) % 2 != 0]

# Prepare the final output format
output = '\n'.join(' '.join(map(str, row)) for row in final_matrix)
print(f"<<<{output}>>>")