# Initial grid setup
grid = [
    [0, 1],
    [0, 0]
]

# Button A controls these positions (0-indexed)
button_a_controls = [(0, 0), (1, 0)]

# Function to toggle the light state
def toggle_light(grid, position):
    x, y = position
    grid[x][y] = 1 - grid[x][y]  # Toggle the light

# Simulate pressing Button A
for position in button_a_controls:
    toggle_light(grid, position)

# Flatten the grid to a single list for the final output
result = [light for row in grid for light in row]
print(result)