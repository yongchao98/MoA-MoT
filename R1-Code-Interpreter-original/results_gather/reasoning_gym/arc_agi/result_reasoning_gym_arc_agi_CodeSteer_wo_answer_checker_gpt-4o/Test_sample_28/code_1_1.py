# Test input grid
input_grid = [
    [6, 6, 0, 6, 6, 6, 6, 6, 6, 6],
    [6, 6, 6, 6, 6, 6, 4, 6, 6, 6],
    [3, 6, 6, 6, 6, 6, 6, 6, 6, 6],
    [6, 6, 6, 6, 1, 6, 6, 6, 6, 6],
    [6, 6, 6, 6, 6, 6, 6, 6, 0, 6],
    [6, 5, 6, 6, 6, 6, 6, 6, 6, 6],
    [6, 6, 6, 6, 6, 6, 6, 6, 3, 6],
    [6, 6, 6, 6, 6, 2, 6, 6, 6, 6],
    [6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
    [6, 6, 6, 6, 6, 6, 6, 9, 6, 6]
]

# Extract non-6 numbers
non_six_numbers = [num for row in input_grid for num in row if num != 6]

# Arrange in spiral order
output_grid = [[0]*3 for _ in range(3)]
positions = [(0, 0), (0, 1), (0, 2), (1, 2), (2, 2), (2, 1), (2, 0), (1, 0), (1, 1)]

for pos, num in zip(positions, non_six_numbers):
    output_grid[pos[0]][pos[1]] = num

# Print the output grid
print("<<<")
for row in output_grid:
    print(" ".join(map(str, row)))
print(">>>")