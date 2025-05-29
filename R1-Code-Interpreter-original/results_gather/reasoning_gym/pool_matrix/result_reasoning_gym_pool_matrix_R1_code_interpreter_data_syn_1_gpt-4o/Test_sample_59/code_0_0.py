import numpy as np

# Define the input matrix
matrix = np.array([
    [2, 2, 4, 4],
    [5, 2, 5, 1],
    [1, 8, 1, 2]
])

# Define the kernel size
kernel_size = 3

# Extract the 3x3 section
section = matrix[:kernel_size, :kernel_size]

# Calculate the average
average = np.mean(section)

# Print the result as a 1x1 matrix
print(f"[[{average:.2f}]]")