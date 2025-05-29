# Initial grid
grid = [
    [0, 0, 0, 1],
    [0, 0, 1, 0],
    [1, 1, 1, 1],
    [1, 0, 1, 0]
]

# Button controls
button_A = [(2, 3), (4, 2), (1, 4), (1, 3), (2, 2), (3, 3), (1, 1), (4, 3)]
button_B = [(1, 4), (3, 2), (4, 1), (2, 2), (3, 3), (1, 2), (3, 4)]
button_C = [(3, 2), (4, 4), (3, 1), (1, 3), (1, 1), (2, 1), (3, 4)]

# Function to toggle lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x-1][y-1] = 1 - grid[x-1][y-1]

# Apply button operations
toggle_lights(grid, button_B)  # Round 1
toggle_lights(grid, button_B)  # Round 2
toggle_lights(grid, button_C)  # Round 3

# Convert grid to a single list
result = [grid[i][j] for i in range(4) for j in range(4)]
print(result)