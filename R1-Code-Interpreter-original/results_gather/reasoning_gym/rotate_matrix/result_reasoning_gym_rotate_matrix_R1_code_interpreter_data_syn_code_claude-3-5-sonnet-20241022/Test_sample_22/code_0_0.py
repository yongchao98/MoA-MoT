import numpy as np

# Input matrix
matrix = np.array([
    [7, 8, 4, 7],
    [7, 6, 3, 7],
    [5, 9, 6, 3],
    [3, 0, 7, 8]
])

# Calculate number of 90-degree rotations needed
degrees = 630
rotations = (degrees // 90) % 4  # We only need to consider rotations modulo 4

# Function to rotate matrix 90 degrees clockwise
def rotate_90_clockwise(mat):
    return np.rot90(mat, k=-1)

# Perform the rotations
result = matrix.copy()
for _ in range(rotations):
    result = rotate_90_clockwise(result)

# Format the output
output = ""
for row in result:
    output += " ".join(map(str, row)) + "\n"
print(output.strip())