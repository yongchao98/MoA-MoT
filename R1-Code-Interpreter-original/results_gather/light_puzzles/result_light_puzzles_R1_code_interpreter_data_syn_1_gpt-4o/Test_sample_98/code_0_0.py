# Initial grid state
grid = [
    [0, 0, 0, 1],
    [0, 0, 1, 0],
    [1, 1, 1, 1],
    [1, 0, 1, 0]
]

# Button control positions
button_A = [(1, 2), (3, 1), (0, 3), (0, 2), (1, 1), (2, 2), (0, 0), (3, 2)]
button_B = [(0, 3), (2, 1), (3, 0), (1, 1), (2, 2), (0, 1), (2, 3)]
button_C = [(2, 1), (3, 3), (2, 0), (0, 2), (0, 0), (1, 0), (2, 3)]

# Function to toggle lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x][y] = 1 - grid[x][y]

# Simulate the rounds
toggle_lights(grid, button_B)  # Round 1: Button B
toggle_lights(grid, button_B)  # Round 2: Button B
toggle_lights(grid, button_C)  # Round 3: Button C

# Convert the grid to a single list
final_state = [grid[i][j] for i in range(4) for j in range(4)]

# Print the final state
print(final_state)