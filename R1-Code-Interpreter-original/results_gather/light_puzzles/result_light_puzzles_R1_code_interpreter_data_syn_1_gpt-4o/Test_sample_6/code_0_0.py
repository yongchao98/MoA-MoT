# Initial grid state
grid = [
    [0, 1, 0, 1, 0],
    [1, 0, 1, 1, 1],
    [0, 0, 0, 0, 0],
    [0, 0, 0, 1, 0],
    [0, 1, 1, 1, 0]
]

# Button control positions
button_A = [(0, 3), (0, 2), (4, 0), (3, 4), (2, 0), (3, 2), (4, 4), (2, 3), (4, 1)]
button_B = [(3, 4), (0, 3), (2, 3), (4, 2), (2, 1), (2, 4), (0, 0), (4, 3), (1, 0), (1, 3), (3, 3)]
button_C = [(0, 3), (2, 3), (0, 2), (2, 0), (4, 4), (0, 1), (2, 2), (4, 0), (3, 0)]

# Function to toggle lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x][y] = 1 - grid[x][y]

# Simulate the rounds
rounds = ['B', 'C', 'B', 'C']
for round in rounds:
    if round == 'A':
        toggle_lights(grid, button_A)
    elif round == 'B':
        toggle_lights(grid, button_B)
    elif round == 'C':
        toggle_lights(grid, button_C)

# Flatten the grid to a single list
final_state = [grid[i][j] for i in range(5) for j in range(5)]
print(final_state)