# Initial grid state
grid = [
    [1, 0, 0, 0, 0],
    [0, 0, 1, 1, 0],
    [1, 0, 0, 0, 1],
    [0, 0, 0, 0, 0],
    [0, 1, 0, 0, 0]
]

# Button control positions (1-based index)
button_A = [(4, 4), (5, 3), (2, 3), (2, 4), (2, 1), (5, 2)]
button_B = [(3, 5), (5, 4), (3, 4), (4, 4)]
button_C = [(2, 1), (2, 5), (4, 2), (1, 5), (3, 4), (5, 5), (5, 2)]

# Convert 1-based index to 0-based index for Python lists
button_A = [(x-1, y-1) for x, y in button_A]
button_B = [(x-1, y-1) for x, y in button_B]
button_C = [(x-1, y-1) for x, y in button_C]

# Function to toggle lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x][y] = 1 - grid[x][y]  # Toggle the light

# Simulate the rounds
toggle_lights(grid, button_A)  # Round 1: Button A
toggle_lights(grid, button_A)  # Round 2: Button A
toggle_lights(grid, button_B)  # Round 3: Button B
toggle_lights(grid, button_C)  # Round 4: Button C

# Flatten the grid to a single list
final_state = [light for row in grid for light in row]

# Print the final state
print(final_state)