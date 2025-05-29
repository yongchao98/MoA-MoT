import numpy as np

# Input matrix
matrix = np.array([
    [6, 5, 6, 8, 2, 8, 8],
    [1, 7, 9, 5, 8, 8, 1]
])

# Add zero padding for the third row to make 3×3 kernels possible
padded = np.pad(matrix, ((0, 1), (0, 0)), mode='constant', constant_values=0)

# Kernel size
k = 3

# Initialize output
cols = matrix.shape[1] // k
output = np.zeros((1, cols))

# Perform average pooling
for j in range(cols):
    # Extract 3×3 region
    region = padded[:, j*k:(j+1)*k]
    # Calculate average
    output[0, j] = np.mean(region)

# Format output to 2 decimal places
formatted_output = np.round(output, 2)

# Print in required format
print('<<<', end='')
print(formatted_output[0][0], formatted_output[0][1], sep=' ', end='')
print('>>>')