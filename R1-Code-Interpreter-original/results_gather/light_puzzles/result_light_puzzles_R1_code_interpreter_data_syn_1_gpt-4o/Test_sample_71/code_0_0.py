# Initial grid configuration
grid = [
    [1, 0, 1],
    [1, 1, 0],
    [0, 0, 1]
]

# Button A controls
button_A = [(2, 1), (0, 2)]

# Button B controls (not used in this scenario)
button_B = [(1, 1), (0, 0), (1, 2), (1, 0)]

# Function to toggle the lights
def toggle_lights(grid, controls):
    for x, y in controls:
        grid[x][y] = 1 - grid[x][y]  # Toggle the light

# Simulate the rounds
# Round 1: Button A
toggle_lights(grid, button_A)

# Round 2: Button A
toggle_lights(grid, button_A)

# Flatten the grid to a single list
result = [grid[i][j] for i in range(3) for j in range(3)]

# Print the result
print(result)