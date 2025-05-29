# Initial grid setup
grid = [
    [1, 0],
    [0, 0]
]

# Button A controls these positions (1-indexed in the problem statement)
button_a_controls = [(1, 1), (1, 0)]

# Function to toggle the light
def toggle_light(grid, x, y):
    grid[x][y] = 1 - grid[x][y]

# Simulate pressing Button A
for x, y in button_a_controls:
    toggle_light(grid, x, y)

# Flatten the grid to a single list for the final output
result = [grid[i][j] for i in range(2) for j in range(2)]

# Print the result
print(result)