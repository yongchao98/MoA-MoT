# Initial grid configuration
grid = [
    [1, 1, 0],
    [1, 1, 1],
    [0, 0, 1]
]

# Button control mappings
button_A = [(1, 1), (1, 2), (0, 2)]
button_B = [(0, 0), (1, 2), (0, 2)]

# Sequence of button presses
rounds = ['B', 'B']

# Function to toggle the light state
def toggle_light(grid, x, y):
    grid[x][y] = 1 - grid[x][y]

# Apply button operations
for round in rounds:
    if round == 'A':
        for x, y in button_A:
            toggle_light(grid, x, y)
    elif round == 'B':
        for x, y in button_B:
            toggle_light(grid, x, y)

# Convert the final grid to a single list
final_state = [grid[i][j] for i in range(3) for j in range(3)]

# Print the final state
print(final_state)