# Initial grid configuration
grid = [
    [1, 0, 1, 0],
    [0, 1, 1, 0],
    [1, 0, 1, 0],
    [1, 1, 0, 1]
]

# Button A controls these positions (1-based index)
button_a_positions = [(4, 4), (3, 4), (2, 2), (1, 3), (2, 4), (1, 4)]

# Button B controls these positions (1-based index)
button_b_positions = [(1, 1), (4, 2), (1, 2), (4, 1)]

# Function to toggle the lights controlled by a button
def toggle_lights(grid, positions):
    for x, y in positions:
        # Convert 1-based index to 0-based index
        grid[x-1][y-1] = 1 - grid[x-1][y-1]

# Simulate the rounds
# Round 1: Button B
toggle_lights(grid, button_b_positions)

# Round 2: Button B
toggle_lights(grid, button_b_positions)

# Round 3: Button A
toggle_lights(grid, button_a_positions)

# Flatten the grid to a single list
final_state = [light for row in grid for light in row]

# Print the final state
print(final_state)