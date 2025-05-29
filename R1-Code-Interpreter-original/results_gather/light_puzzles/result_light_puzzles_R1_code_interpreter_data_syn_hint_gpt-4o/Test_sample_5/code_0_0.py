# Initial grid state
grid = [
    [1, 1, 0, 0, 1],
    [0, 1, 0, 1, 0],
    [0, 1, 1, 0, 0],
    [1, 0, 0, 0, 0],
    [1, 1, 0, 1, 1]
]

# Button control positions
button_A = [(3, 4), (2, 5), (2, 4), (5, 1), (1, 5), (5, 4), (3, 1), (1, 4), (4, 5), (5, 2)]
button_B = [(4, 5), (1, 3), (5, 1), (5, 5)]
button_C = [(4, 3), (5, 4), (4, 1), (5, 1), (3, 4), (2, 1), (5, 3), (2, 5), (5, 2)]

# Function to toggle lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x-1][y-1] = 1 - grid[x-1][y-1]

# Simulate rounds
toggle_lights(grid, button_A)  # Round 1: Button A
toggle_lights(grid, button_A)  # Round 2: Button A
toggle_lights(grid, button_C)  # Round 3: Button C
toggle_lights(grid, button_C)  # Round 4: Button C

# Flatten the grid to a single list
final_state = [grid[i][j] for i in range(5) for j in range(5)]

# Print the final state
print(final_state)