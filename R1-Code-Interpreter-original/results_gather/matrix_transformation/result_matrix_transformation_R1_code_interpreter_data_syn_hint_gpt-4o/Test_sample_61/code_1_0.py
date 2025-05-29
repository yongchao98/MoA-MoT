# Rotated matrix in 2D form
rotated_matrix_2d = [
    ['9', 'I', 'G'],
    ['8', 'O', 'P'],
    ['X', 'E', '7'],
    ['X', '1', 'L']
]

# Flatten the matrix and join elements with commas
formatted_result = ','.join([item for sublist in rotated_matrix_2d for item in sublist])

# Print the formatted result
print(f"<<<{formatted_result}>>>")