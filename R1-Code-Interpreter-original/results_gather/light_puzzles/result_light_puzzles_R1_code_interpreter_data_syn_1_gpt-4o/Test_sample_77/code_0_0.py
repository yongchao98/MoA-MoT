# Initial grid configuration
grid = [
    [0, 1, 1],
    [0, 0, 0],
    [0, 0, 1]
]

# Button A controls
button_A = [(2, 0), (2, 1)]  # 0-indexed positions

# Button B controls
button_B = [(2, 1), (2, 0)]  # 0-indexed positions

# Function to toggle lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x][y] = 1 - grid[x][y]  # Toggle the light

# Round 1: Button A
toggle_lights(grid, button_A)

# Round 2: Button B
toggle_lights(grid, button_B)

# Flatten the grid to a single list
result = [grid[i][j] for i in range(3) for j in range(3)]

# Print the result
print(result)