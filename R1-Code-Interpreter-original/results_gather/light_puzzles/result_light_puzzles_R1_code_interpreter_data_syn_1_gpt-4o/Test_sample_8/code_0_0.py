# Initial grid state
grid = [
    [0, 0, 1, 1, 1],
    [1, 1, 0, 0, 1],
    [0, 1, 1, 1, 1],
    [1, 0, 0, 0, 1],
    [0, 1, 0, 0, 1]
]

# Button control positions (1-based index)
button_A = [(2, 4), (1, 3), (4, 2), (2, 1), (4, 3), (5, 1), (2, 3)]
button_B = [(1, 2), (3, 4), (2, 1), (1, 3), (2, 4)]
button_C = [(4, 3), (3, 4), (1, 3), (4, 4), (2, 1)]

# Convert to 0-based index for easier access in Python
button_A = [(x-1, y-1) for x, y in button_A]
button_B = [(x-1, y-1) for x, y in button_B]
button_C = [(x-1, y-1) for x, y in button_C]

# Sequence of button presses
rounds = ['C', 'C', 'A', 'C']

# Function to toggle lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x][y] = 1 - grid[x][y]  # Toggle between 0 and 1

# Apply each round
for round in rounds:
    if round == 'A':
        toggle_lights(grid, button_A)
    elif round == 'B':
        toggle_lights(grid, button_B)
    elif round == 'C':
        toggle_lights(grid, button_C)

# Flatten the grid to a single list
final_state = [light for row in grid for light in row]

# Print the final state
print(final_state)