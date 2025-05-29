# Initial grid configuration
grid = [
    [1, 0, 0, 0],
    [1, 1, 0, 1],
    [0, 1, 1, 1],
    [0, 1, 0, 1]
]

# Button control positions
button_A_positions = [(2, 3), (1, 2), (3, 2), (2, 2), (2, 4), (4, 3), (1, 3)]
button_B_positions = [(2, 3), (4, 4), (2, 4), (3, 1), (4, 2), (4, 3), (2, 2), (2, 1)]

# Function to toggle lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x-1][y-1] = 1 - grid[x-1][y-1]  # Toggle the light

# Simulate the rounds
toggle_lights(grid, button_B_positions)  # Round 1: Button B
toggle_lights(grid, button_A_positions)  # Round 2: Button A
toggle_lights(grid, button_A_positions)  # Round 3: Button A

# Convert the grid to a single list of integers
final_state = [grid[i][j] for i in range(4) for j in range(4)]
print(final_state)