# Initial grid configuration
grid = [
    [1, 0, 0],
    [1, 1, 1],
    [0, 1, 1]
]

# Button A controls
button_A = [(1, 2), (2, 0), (1, 1)]

# Button B controls
button_B = [(2, 1), (0, 0), (2, 0)]

# Function to toggle lights
def toggle_lights(grid, controls):
    for x, y in controls:
        grid[x][y] = 1 - grid[x][y]  # Toggle the light

# Round 1: Button A
toggle_lights(grid, button_A)

# Round 2: Button B
toggle_lights(grid, button_B)

# Flatten the grid to a single list
final_state = [grid[i][j] for i in range(3) for j in range(3)]

# Print the final state
print(final_state)