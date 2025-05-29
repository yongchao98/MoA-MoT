# Initial grid
grid = [
    [1, 0],
    [0, 1]
]

# Button A controls these positions
button_a_positions = [(0, 0), (1, 0)]

# Function to toggle the light
def toggle_light(grid, position):
    x, y = position
    grid[x][y] = 1 - grid[x][y]  # Toggle between 0 and 1

# Simulate pressing Button A
for position in button_a_positions:
    toggle_light(grid, position)

# Flatten the grid to a single list for the final output
result = [grid[i][j] for i in range(2) for j in range(2)]
print(result)