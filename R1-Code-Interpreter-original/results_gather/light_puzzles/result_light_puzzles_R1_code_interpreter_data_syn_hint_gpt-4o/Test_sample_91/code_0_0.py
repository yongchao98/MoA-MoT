# Initial grid setup
grid = [
    [0, 1, 0, 0],
    [1, 0, 0, 1],
    [1, 1, 1, 1],
    [1, 1, 1, 0]
]

# Button control mappings
button_A_controls = [(0, 3), (0, 2), (0, 1), (2, 1)]
button_B_controls = [(0, 3), (2, 0), (1, 2), (3, 0), (1, 0)]

# Function to toggle lights
def toggle_lights(grid, controls):
    for x, y in controls:
        grid[x][y] = 1 - grid[x][y]

# Round 1: Press Button B
toggle_lights(grid, button_B_controls)

# Round 2: Press Button A
toggle_lights(grid, button_A_controls)

# Round 3: Press Button A
toggle_lights(grid, button_A_controls)

# Flatten the grid to a single list
result = [grid[i][j] for i in range(4) for j in range(4)]
print(result)