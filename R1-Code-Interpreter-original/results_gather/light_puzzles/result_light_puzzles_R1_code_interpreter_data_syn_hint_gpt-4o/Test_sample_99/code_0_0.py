# Initial grid setup
grid = [
    [1, 0, 1, 0],
    [1, 1, 1, 0],
    [0, 1, 0, 0],
    [1, 1, 1, 0]
]

# Button control positions
button_A = [(0, 3), (0, 0), (1, 1), (0, 2), (3, 1)]
button_B = [(3, 3), (0, 0), (2, 1), (3, 2), (0, 3), (3, 1)]
button_C = [(0, 3), (1, 0), (0, 2), (1, 3), (3, 2), (3, 0)]

# Function to toggle lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x][y] = 1 - grid[x][y]  # Toggle the light

# Simulate the rounds
toggle_lights(grid, button_A)  # Round 1: Button A
toggle_lights(grid, button_A)  # Round 2: Button A
toggle_lights(grid, button_C)  # Round 3: Button C

# Flatten the grid to a single list
final_state = [grid[i][j] for i in range(4) for j in range(4)]

# Print the final state
print(final_state)