# Initial grid state
grid = [
    [1, 0, 0, 1],
    [1, 1, 1, 1],
    [1, 1, 1, 0],
    [1, 0, 0, 1]
]

# Button control positions
button_A = [(4, 3), (4, 4), (1, 1), (1, 4)]
button_B = [(4, 1), (2, 1), (2, 4), (2, 3), (3, 4), (4, 4)]
button_C = [(4, 4), (3, 2), (1, 3), (3, 4), (3, 3), (4, 2), (2, 1), (1, 4)]

# Convert 1-based positions to 0-based for Python indexing
button_A = [(x-1, y-1) for x, y in button_A]
button_B = [(x-1, y-1) for x, y in button_B]
button_C = [(x-1, y-1) for x, y in button_C]

# Function to toggle lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x][y] = 1 - grid[x][y]  # Toggle between 0 and 1

# Simulate the rounds
toggle_lights(grid, button_B)  # Round 1: Button B
toggle_lights(grid, button_A)  # Round 2: Button A
toggle_lights(grid, button_C)  # Round 3: Button C

# Flatten the grid to a single list
final_state = [grid[i][j] for i in range(4) for j in range(4)]

# Print the final state
print(final_state)