# Initial grid configuration
grid = [
    [0, 0, 1, 1, 1],
    [1, 1, 0, 0, 1],
    [0, 1, 1, 1, 1],
    [1, 0, 0, 0, 1],
    [0, 1, 0, 0, 1]
]

# Button controls
button_controls = {
    'A': [(2, 4), (1, 3), (4, 2), (2, 1), (4, 3), (5, 1), (2, 3)],
    'B': [(1, 2), (3, 4), (2, 1), (1, 3), (2, 4)],
    'C': [(4, 3), (3, 4), (1, 3), (4, 4), (2, 1)]
}

# Rounds of button presses
rounds = ['C', 'C', 'A', 'C']

# Function to toggle the light at a given position
def toggle_light(grid, x, y):
    grid[x][y] = 1 - grid[x][y]

# Simulate the rounds
for button in rounds:
    for x, y in button_controls[button]:
        toggle_light(grid, x - 1, y - 1)  # Convert to 0-based index

# Flatten the grid to a single list
final_state = [light for row in grid for light in row]

# Print the final state
print(final_state)