# Input matrix
matrix = [
    [1, 9],
    [8, 5],
    [0, 1],
    [4, 1],
    [9, 5],
    [9, 8]
]

kernel_size = 2
rows = len(matrix)
cols = len(matrix[0])
output = []

# Process the matrix in 2×2 blocks
for i in range(0, rows, kernel_size):
    sum_elements = 0
    # Sum all elements in the current 2×2 block
    for r in range(i, min(i + kernel_size, rows)):
        for c in range(cols):
            sum_elements += matrix[r][c]
    # Calculate average and append to output
    average = sum_elements / (kernel_size * cols)
    output.append(round(average, 2))

# Print the result in the required format
print("<<<")
for val in output:
    print(f"{val:.2f}")
print(">>>")