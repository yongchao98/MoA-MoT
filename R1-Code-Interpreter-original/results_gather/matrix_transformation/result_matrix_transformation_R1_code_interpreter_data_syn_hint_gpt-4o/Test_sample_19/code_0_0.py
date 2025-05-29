# Original matrix
matrix = [
    ["W", "T"],
    ["Z", "B"]
]

# Rotating the matrix 90 degrees clockwise
rotated_matrix = list(zip(*matrix[::-1]))

# Flatten the rotated matrix to match the required output format
flattened_result = [item for sublist in rotated_matrix for item in sublist]

# Print the result in the required format
print(flattened_result)