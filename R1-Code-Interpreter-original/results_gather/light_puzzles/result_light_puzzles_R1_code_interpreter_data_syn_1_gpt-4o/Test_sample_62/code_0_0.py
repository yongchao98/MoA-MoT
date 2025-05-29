# Initial grid setup
grid = [
    [1, 0],
    [0, 0]
]

# Button A affects these positions (1-indexed): ["2", "2"], ["2", "1"]
# Convert to 0-indexed for Python: (1, 1), (1, 0)
positions_to_toggle = [(1, 1), (1, 0)]

# Function to toggle the light
def toggle_light(grid, position):
    x, y = position
    grid[x][y] = 1 - grid[x][y]  # Toggle: 1 becomes 0, 0 becomes 1

# Simulate pressing Button A
for position in positions_to_toggle:
    toggle_light(grid, position)

# Flatten the grid to a single list for the final output
result = [light for row in grid for light in row]

# Print the result
print(result)