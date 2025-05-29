# Initial grid configuration
grid = [
    [0, 0, 1],
    [0, 0, 0],
    [0, 1, 1]
]

# Function to toggle the state of a light
def toggle_light(grid, x, y):
    grid[x][y] = 1 - grid[x][y]

# Button A controls
button_a_controls = [(0, 2), (1, 1), (0, 1), (1, 2)]

# Button B controls
button_b_controls = [(2, 2), (0, 1), (1, 1)]

# Simulate button presses
def press_button(grid, controls):
    for x, y in controls:
        toggle_light(grid, x, y)

# Round 1: Button A
press_button(grid, button_a_controls)

# Round 2: Button A
press_button(grid, button_a_controls)

# Flatten the grid to a single list
final_state = [grid[i][j] for i in range(3) for j in range(3)]

# Print the final state
print(final_state)