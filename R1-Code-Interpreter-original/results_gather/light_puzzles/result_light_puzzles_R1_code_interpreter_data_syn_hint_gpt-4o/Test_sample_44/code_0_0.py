# Initial grid state
grid = [
    [1, 1, 0, 1, 1],
    [1, 0, 0, 1, 1],
    [1, 1, 0, 1, 1],
    [0, 0, 0, 0, 1],
    [0, 0, 0, 0, 1]
]

# Button control positions
button_A_positions = [(5, 3), (3, 5), (4, 4), (1, 4), (2, 3), (2, 1), (4, 1), (1, 5), (5, 4), (3, 1), (5, 5), (1, 1)]
button_B_positions = [(5, 5), (3, 3), (1, 1), (3, 4), (1, 4), (2, 1), (3, 1), (2, 2)]
button_C_positions = [(3, 5), (1, 5), (4, 2), (5, 3)]

# Convert positions to zero-indexed
button_A_positions = [(x-1, y-1) for x, y in button_A_positions]
button_B_positions = [(x-1, y-1) for x, y in button_B_positions]
button_C_positions = [(x-1, y-1) for x, y in button_C_positions]

# Function to toggle lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x][y] = 1 - grid[x][y]

# Simulate rounds
toggle_lights(grid, button_A_positions)  # Round 1: Button A
toggle_lights(grid, button_B_positions)  # Round 2: Button B
toggle_lights(grid, button_B_positions)  # Round 3: Button B
toggle_lights(grid, button_A_positions)  # Round 4: Button A

# Flatten the grid to a single list
final_state = [grid[i][j] for i in range(5) for j in range(5)]

# Print the final state
print(final_state)