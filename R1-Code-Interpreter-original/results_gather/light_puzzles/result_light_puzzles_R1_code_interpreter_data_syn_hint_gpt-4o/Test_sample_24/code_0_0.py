# Initial grid setup
grid = [
    [1, 1, 1, 0, 1],
    [0, 0, 1, 0, 1],
    [1, 0, 0, 0, 0],
    [0, 0, 0, 0, 0],
    [1, 0, 0, 1, 0]
]

# Button B controls these positions (1-based index)
button_b_positions = [
    (5, 5), (3, 4), (5, 3), (5, 4), (5, 1), (3, 5), (1, 5), (2, 5)
]

# Convert to 0-based index for Python
button_b_positions = [(x-1, y-1) for x, y in button_b_positions]

# Function to toggle the light at a given position
def toggle_light(grid, x, y):
    grid[x][y] = 1 - grid[x][y]

# Simulate pressing Button B four times
for _ in range(4):
    for x, y in button_b_positions:
        toggle_light(grid, x, y)

# Flatten the grid to a single list
final_state = [grid[i][j] for i in range(5) for j in range(5)]

# Print the final state
print(final_state)