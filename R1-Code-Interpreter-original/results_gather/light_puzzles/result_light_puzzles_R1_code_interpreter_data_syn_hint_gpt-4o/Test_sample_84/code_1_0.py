# Initial grid setup
grid = [
    [0, 1, 0, 0],
    [1, 1, 1, 0],
    [0, 1, 1, 0],
    [0, 1, 0, 0]
]

# Button B controls these positions (1-based index)
button_b_positions = [(3, 3), (3, 2), (1, 4), (3, 1)]

# Function to toggle the light at a given position
def toggle_light(grid, x, y):
    grid[x][y] = 1 - grid[x][y]

# Simulate the button presses
for _ in range(3):  # Three rounds of pressing Button B
    for x, y in button_b_positions:
        toggle_light(grid, x-1, y-1)  # Convert to 0-based index

# Flatten the grid to a single list
result = [grid[i][j] for i in range(4) for j in range(4)]

# Print the result
print(result)