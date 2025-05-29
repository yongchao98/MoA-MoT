# Initial grid state
grid = [
    [1, 0, 0, 0, 0],
    [1, 1, 0, 1, 1],
    [0, 0, 0, 1, 1],
    [1, 1, 1, 0, 0],
    [1, 0, 1, 0, 0]
]

# Button control positions
button_A = [(2, 4), (4, 2), (2, 3), (3, 5), (4, 4), (1, 3), (1, 1), (4, 3), (3, 1), (5, 4), (4, 1)]
button_B = [(1, 2), (3, 5), (4, 2), (5, 3), (3, 4), (3, 1), (3, 3)]
button_C = [(5, 5), (5, 1), (4, 2), (4, 1), (3, 1), (2, 3), (3, 4), (5, 2), (1, 5), (1, 1), (3, 5)]

# Convert to zero-indexed
button_A = [(x-1, y-1) for x, y in button_A]
button_B = [(x-1, y-1) for x, y in button_B]
button_C = [(x-1, y-1) for x, y in button_C]

# Function to toggle lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x][y] = 1 - grid[x][y]

# Simulate rounds
rounds = ['C', 'A', 'C', 'C']
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