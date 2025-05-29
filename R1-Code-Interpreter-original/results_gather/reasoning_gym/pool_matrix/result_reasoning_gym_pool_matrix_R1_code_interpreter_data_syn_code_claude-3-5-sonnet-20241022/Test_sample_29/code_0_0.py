import numpy as np

# Input matrix
matrix = [
    [3, 2, 2, 0],
    [2, 1, 7, 5],
    [3, 7, 0, 9],
    [3, 4, 6, 2],
    [3, 5, 8, 2],
    [6, 3, 2, 9],
    [1, 1, 6, 5],
    [0, 2, 6, 2],
    [2, 9, 4, 5]
]

kernel_size = 3

# Convert to numpy array for easier manipulation
matrix = np.array(matrix)
rows, cols = matrix.shape

# Calculate output dimensions
out_rows = rows // kernel_size
out_cols = cols // kernel_size

# Initialize output matrix
output = np.zeros((out_rows, out_cols))

# Perform average pooling
for i in range(out_rows):
    for j in range(out_cols):
        # Extract the current kernel region
        region = matrix[i*kernel_size:(i+1)*kernel_size, 
                       j*kernel_size:(j+1)*kernel_size]
        # Calculate average
        output[i,j] = np.mean(region)

# Format output to 2 decimal places
formatted_output = np.round(output, 2)

# Print result in required format
print("<<<")
for row in formatted_output:
    print(" ".join(f"{x:.2f}" for x in row))
print(">>>")