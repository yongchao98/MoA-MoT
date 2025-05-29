# Initial grid
grid = [
    [1, 0],
    [0, 0]
]

# Button A controls these positions
button_a_controls = [(1, 0), (1, 1)]  # (row, column) in 0-indexed

# Toggle the lights controlled by Button A
for row, col in button_a_controls:
    grid[row][col] = 1 - grid[row][col]  # Toggle the light

# Flatten the grid to a single list
result = [light for row in grid for light in row]

# Print the result
print(result)