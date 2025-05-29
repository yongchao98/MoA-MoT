# Original matrix
matrix = [
    ["T", "S", "I", "X"],
    ["5", "5", "3", "E"],
    ["I", "S", "7", "H"]
]

# Rotating the matrix 90 degrees clockwise
rotated_matrix = list(zip(*matrix[::-1]))

# Formatting the output
result = [",".join(row) for row in rotated_matrix]
print("<<<" + ",".join(result) + ">>>")