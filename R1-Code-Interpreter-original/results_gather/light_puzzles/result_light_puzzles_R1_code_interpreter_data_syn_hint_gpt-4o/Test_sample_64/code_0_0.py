# Initial grid configuration
grid = [
    [1, 1, 1],
    [1, 0, 1],
    [0, 1, 0]
]

# Button controls
button_A_controls = [(2, 2), (0, 0)]
button_B_controls = [(2, 0), (0, 2)]

# Function to toggle the light at a given position
def toggle_light(grid, x, y):
    grid[x][y] = 1 - grid[x][y]

# Simulate the rounds
rounds = ['A', 'A']

for round in rounds:
    if round == 'A':
        for x, y in button_A_controls:
            toggle_light(grid, x, y)
    elif round == 'B':
        for x, y in button_B_controls:
            toggle_light(grid, x, y)

# Flatten the grid to a single list
final_state = [grid[i][j] for i in range(3) for j in range(3)]

# Print the final state
print(final_state)