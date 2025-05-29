# Initial grid configuration
grid = [
    [1, 1, 1],
    [1, 0, 1],
    [0, 1, 0]
]

# Button A controls
button_A = [(2, 2), (0, 0)]

# Button B controls
button_B = [(2, 0), (0, 2)]

# Function to toggle the lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x][y] = 1 - grid[x][y]  # Toggle the light

# Simulate the rounds
# Round 1: Button A
toggle_lights(grid, button_A)

# Round 2: Button A
toggle_lights(grid, button_A)

# Flatten the grid to a single list
final_state = [grid[i][j] for i in range(3) for j in range(3)]

# Print the final state
print(final_state)