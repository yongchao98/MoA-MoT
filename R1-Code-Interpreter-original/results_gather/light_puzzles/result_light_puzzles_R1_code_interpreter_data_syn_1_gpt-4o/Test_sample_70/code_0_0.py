# Initial grid configuration
grid = [
    [0, 0, 0],
    [0, 1, 0],
    [1, 0, 1]
]

# Button A controls
button_A = [(1, 0), (1, 1), (2, 1)]

# Button B controls
button_B = [(1, 1), (2, 0), (0, 1)]

# Function to toggle the lights
def toggle_lights(grid, controls):
    for x, y in controls:
        grid[x][y] = 1 - grid[x][y]  # Toggle the light

# Simulate the rounds
# Round 1: Button B
toggle_lights(grid, button_B)

# Round 2: Button B
toggle_lights(grid, button_B)

# Convert the grid to a single list
result = [grid[i][j] for i in range(3) for j in range(3)]

# Print the result
print(result)