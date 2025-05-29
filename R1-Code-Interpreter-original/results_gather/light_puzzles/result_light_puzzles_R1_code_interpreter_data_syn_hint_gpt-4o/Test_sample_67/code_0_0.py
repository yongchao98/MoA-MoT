# Initial grid setup
n = 3
grid = [
    [0, 0, 1],
    [0, 0, 0],
    [0, 1, 1]
]

# Button control positions (1-based index)
button_A_positions = [(1, 3), (2, 2), (1, 2), (2, 3)]
button_B_positions = [(3, 3), (1, 2), (2, 2)]

# Convert to 0-based index for easier manipulation
button_A_positions = [(x-1, y-1) for x, y in button_A_positions]
button_B_positions = [(x-1, y-1) for x, y in button_B_positions]

# Function to toggle lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x][y] = 1 - grid[x][y]  # Toggle between 0 and 1

# Simulate the rounds
toggle_lights(grid, button_A_positions)  # Round 1: Button A
toggle_lights(grid, button_A_positions)  # Round 2: Button A

# Flatten the grid to a single list
final_state = [grid[i][j] for i in range(n) for j in range(n)]

# Print the final state
print(final_state)