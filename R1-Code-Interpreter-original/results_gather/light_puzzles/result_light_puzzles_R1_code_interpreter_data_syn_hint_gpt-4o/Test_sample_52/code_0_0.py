# Initial grid
grid = [
    [1, 0],
    [1, 1]
]

# Button A controls these positions (0-based index)
button_a_controls = [(1, 0), (0, 0)]

# Function to toggle the light
def toggle_light(grid, position):
    x, y = position
    grid[x][y] = 1 - grid[x][y]  # Toggle the light

# Simulate pressing Button A
for position in button_a_controls:
    toggle_light(grid, position)

# Flatten the grid to a single list for the final output
result = [grid[i][j] for i in range(2) for j in range(2)]
print(result)