# Initial grid state
grid = [
    [1, 1, 1, 1, 1],
    [0, 1, 1, 0, 0],
    [0, 1, 1, 1, 0],
    [1, 0, 0, 0, 0],
    [1, 1, 0, 1, 0]
]

# Button control positions (1-based index)
button_A = [(1, 5), (2, 5), (5, 3), (4, 3), (2, 3), (4, 4), (1, 4)]
button_B = [(2, 4), (3, 3), (2, 3)]
button_C = [(1, 3), (3, 5)]

# Convert to 0-based index for easier manipulation
button_A = [(x-1, y-1) for x, y in button_A]
button_B = [(x-1, y-1) for x, y in button_B]
button_C = [(x-1, y-1) for x, y in button_C]

# Function to toggle lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x][y] = 1 - grid[x][y]  # Toggle between 0 and 1

# Simulate each round
toggle_lights(grid, button_A)  # Round 1: Button A
toggle_lights(grid, button_A)  # Round 2: Button A
toggle_lights(grid, button_C)  # Round 3: Button C
toggle_lights(grid, button_C)  # Round 4: Button C

# Flatten the grid to a single list
final_state = [cell for row in grid for cell in row]

# Print the final state
print(final_state)