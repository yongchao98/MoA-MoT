# Initial grid setup
grid = [
    [1, 0],
    [1, 0]
]

# Button A controls these positions
button_a_positions = [(1, 0), (1, 1)]

# Function to toggle the light
def toggle_light(grid, position):
    x, y = position
    grid[x][y] = 1 - grid[x][y]  # Toggle the light

# Simulate pressing Button A
for position in button_a_positions:
    toggle_light(grid, position)

# Flatten the grid to a single list for the final output
result = [light for row in grid for light in row]

# Print the result
print(result)