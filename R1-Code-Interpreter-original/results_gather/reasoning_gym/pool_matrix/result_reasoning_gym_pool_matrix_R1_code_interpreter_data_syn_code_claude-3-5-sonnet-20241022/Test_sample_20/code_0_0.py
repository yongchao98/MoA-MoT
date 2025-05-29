import numpy as np

# Input matrix
matrix = np.array([
    [5, 7],
    [6, 7],
    [0, 5],
    [3, 3],
    [1, 3],
    [1, 5]
])

kernel_size = 3
rows, cols = matrix.shape
output_rows = rows // kernel_size

# Initialize output matrix
output = np.zeros((output_rows, cols))

# Perform max pooling
for i in range(output_rows):
    start_row = i * kernel_size
    end_row = start_row + kernel_size
    region = matrix[start_row:end_row, :]
    output[i] = np.max(region, axis=0)

# Format output to 2 decimal places
formatted_output = np.round(output, 2)

# Print result in required format
print("<<<")
for row in formatted_output:
    print(" ".join(f"{x:.2f}" for x in row))
print(">>>")