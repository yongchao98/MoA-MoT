# Initial grid configuration
grid = [
    [1, 0, 1, 1, 1],
    [0, 1, 0, 0, 1],
    [0, 1, 1, 0, 0],
    [1, 0, 1, 1, 0],
    [1, 1, 1, 0, 1]
]

# Button control positions
button_A = [(3, 5), (2, 3), (2, 2), (4, 5), (4, 2), (1, 4), (3, 4), (1, 1), (3, 3), (2, 5)]
button_B = [(5, 1), (3, 4), (2, 5), (1, 1), (5, 4), (2, 2), (4, 5), (3, 5)]
button_C = [(2, 1), (5, 2), (1, 1), (4, 1), (2, 5), (3, 3), (2, 2), (3, 2), (5, 4), (5, 3), (4, 2), (3, 5)]

# Function to toggle lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x-1][y-1] = 1 - grid[x-1][y-1]

# Simulate the rounds
toggle_lights(grid, button_C)  # Round 1: Button C
toggle_lights(grid, button_C)  # Round 2: Button C
toggle_lights(grid, button_A)  # Round 3: Button A
toggle_lights(grid, button_B)  # Round 4: Button B

# Flatten the grid to a single list
final_state = [grid[i][j] for i in range(5) for j in range(5)]

# Print the final state
print(final_state)