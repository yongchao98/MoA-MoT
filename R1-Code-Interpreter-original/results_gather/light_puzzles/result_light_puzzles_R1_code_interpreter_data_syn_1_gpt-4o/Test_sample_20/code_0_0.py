# Initial grid state
grid = [
    [1, 1, 1, 1, 1],
    [0, 1, 1, 0, 0],
    [0, 1, 1, 1, 0],
    [1, 0, 0, 0, 0],
    [1, 1, 0, 1, 0]
]

# Button control mappings
button_A = [(0, 4), (1, 4), (4, 2), (3, 2), (1, 2), (3, 3), (0, 3)]
button_B = [(1, 3), (2, 2), (1, 2)]
button_C = [(0, 2), (2, 4)]

# Function to toggle lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x][y] = 1 - grid[x][y]

# Simulate the rounds
toggle_lights(grid, button_A)  # Round 1: Button A
toggle_lights(grid, button_A)  # Round 2: Button A
toggle_lights(grid, button_C)  # Round 3: Button C
toggle_lights(grid, button_C)  # Round 4: Button C

# Flatten the grid to a single list
final_state = [grid[i][j] for i in range(5) for j in range(5)]

# Print the final state
print(final_state)