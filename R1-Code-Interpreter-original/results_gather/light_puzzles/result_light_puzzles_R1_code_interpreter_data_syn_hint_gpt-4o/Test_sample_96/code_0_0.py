# Initial grid setup
grid = [
    [1, 0, 1, 1],
    [0, 1, 0, 1],
    [1, 0, 0, 1],
    [0, 0, 1, 1]
]

# Button control mappings
button_A_controls = [(3, 1), (0, 1), (2, 1), (2, 2), (0, 0), (2, 0)]
button_B_controls = [(2, 2), (0, 0), (3, 3), (1, 0)]
button_C_controls = [(2, 1), (2, 2), (2, 3), (0, 3), (0, 2), (3, 1)]

# Function to toggle lights
def toggle_lights(grid, controls):
    for x, y in controls:
        grid[x][y] = 1 - grid[x][y]  # Toggle the light

# Simulate the rounds
for _ in range(3):  # Three rounds of pressing Button A
    toggle_lights(grid, button_A_controls)

# Flatten the grid to a single list
final_state = [grid[i][j] for i in range(4) for j in range(4)]

# Print the final state
print(final_state)