import numpy as np

def apply_kernel(grid, kernel, threshold):
    rows, cols = len(grid), len(grid[0])
    output_grid = [row[:] for row in grid]  # Copy the input grid
    kernel_size = len(kernel)
    offset = kernel_size // 2

    for i in range(offset, rows - offset):
        for j in range(offset, cols - offset):
            # Extract the subgrid
            subgrid = [grid[i + ki - offset][j + kj - offset] for ki in range(kernel_size) for kj in range(kernel_size)]
            # Apply the kernel
            result = sum(k * s for k, s in zip(kernel.flatten(), subgrid))
            # Check if the pattern is detected
            if result >= threshold:
                output_grid[i][j] = 5

    return output_grid

# Define the kernel for detecting the cross pattern
kernel = np.array([
    [0, 1, 0],
    [1, 0, 1],
    [0, 1, 0]
])

# Test input grid
input_grid = [
    [4, 4, 4, 4, 4, 4, 8, 8, 8, 8, 8, 4, 4, 4, 4, 4],
    [8, 8, 8, 8, 4, 4, 8, 4, 4, 8, 4, 4, 4, 4, 4, 4],
    [4, 4, 4, 8, 4, 4, 4, 4, 4, 8, 4, 4, 4, 4, 4, 4],
    [4, 4, 4, 4, 4, 4, 4, 4, 8, 8, 8, 8, 8, 8, 8, 4],
    [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4],
    [8, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4],
    [8, 8, 8, 8, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4],
    [8, 4, 4, 8, 4, 4, 4, 4, 8, 8, 8, 8, 8, 8, 4, 4],
    [8, 4, 4, 8, 4, 4, 4, 4, 4, 4, 8, 4, 4, 8, 4, 4],
    [8, 4, 4, 8, 4, 4, 4, 4, 4, 4, 8, 4, 4, 8, 4, 4],
    [8, 4, 4, 8, 4, 4, 4, 4, 4, 4, 8, 8, 8, 8, 4, 4],
    [8, 8, 8, 8, 4, 4, 8, 8, 8, 4, 4, 4, 4, 8, 4, 4],
    [4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 4, 4, 4, 4, 4, 4],
    [4, 4, 4, 4, 4, 8, 4, 4, 8, 4, 4, 4, 4, 4, 4, 8],
    [4, 4, 4, 4, 4, 8, 8, 8, 8, 4, 4, 4, 4, 4, 4, 8]
]

# Convert input grid to numpy array for easier manipulation
input_grid_np = np.array(input_grid)

# Apply the kernel with a threshold of 3 (since we expect 3 '8's around a '4')
output_grid = apply_kernel(input_grid_np, kernel, threshold=3)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))