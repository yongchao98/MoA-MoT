# Initial grid state
grid = [
    [1, 1, 1, 0, 0],
    [1, 0, 1, 1, 1],
    [0, 0, 1, 0, 0],
    [1, 0, 0, 0, 1],
    [0, 1, 0, 0, 1]
]

# Button control positions
button_A = [(3, 4), (3, 0), (1, 4), (4, 3)]
button_B = [(1, 1), (3, 3), (0, 1), (4, 3)]
button_C = [(3, 0), (3, 4), (4, 3), (1, 1), (4, 4), (0, 1), (0, 2)]

# Rounds of button presses
rounds = ['C', 'C', 'B', 'B']

# Function to toggle the lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x][y] = 1 - grid[x][y]

# Simulate each round
for round in rounds:
    if round == 'A':
        toggle_lights(grid, button_A)
    elif round == 'B':
        toggle_lights(grid, button_B)
    elif round == 'C':
        toggle_lights(grid, button_C)

# Flatten the grid to a single list
final_state = [grid[i][j] for i in range(5) for j in range(5)]

# Print the final state
print(final_state)