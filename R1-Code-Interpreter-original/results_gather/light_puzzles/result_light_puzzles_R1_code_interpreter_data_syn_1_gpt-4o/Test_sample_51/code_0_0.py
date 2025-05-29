# Initial grid
grid = [
    [1, 1],
    [1, 0]
]

# Positions controlled by Button A
button_a_positions = [(1, 1), (0, 1)]

# Function to toggle the light state
def toggle_light(grid, position):
    x, y = position
    grid[x][y] = 1 - grid[x][y]  # Toggle between 0 and 1

# Simulate pressing Button A
for position in button_a_positions:
    toggle_light(grid, position)

# Flatten the grid to a single list
result = [grid[i][j] for i in range(2) for j in range(2)]

# Print the result
print(result)