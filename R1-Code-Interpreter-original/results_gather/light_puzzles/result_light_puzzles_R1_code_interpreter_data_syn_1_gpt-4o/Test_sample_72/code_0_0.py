# Initial grid configuration
grid = [
    [0, 1, 0],
    [1, 0, 1],
    [0, 1, 0]
]

# Button configurations
button_A = [(1, 2), (1, 0)]
button_B = [(1, 0), (1, 1), (0, 0)]

# Function to toggle the lights controlled by a button
def toggle_lights(grid, button):
    for x, y in button:
        grid[x][y] = 1 - grid[x][y]  # Toggle the light

# Simulate the rounds
toggle_lights(grid, button_A)  # Round 1: Button A
toggle_lights(grid, button_B)  # Round 2: Button B

# Flatten the grid to a single list
final_state = [grid[i][j] for i in range(3) for j in range(3)]

# Print the final state
print(final_state)