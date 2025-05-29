# Initial grid setup
grid = [
    [1, 0, 1, 1],
    [0, 1, 0, 1],
    [1, 0, 0, 1],
    [0, 0, 1, 1]
]

# Button control positions (1-based index)
button_A_positions = [(4, 2), (1, 2), (3, 2), (3, 3), (1, 1), (3, 1)]
button_B_positions = [(3, 3), (1, 1), (4, 4), (2, 1)]
button_C_positions = [(3, 2), (3, 3), (3, 4), (1, 4), (1, 3), (4, 2)]

# Function to toggle lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x-1][y-1] = 1 - grid[x-1][y-1]  # Toggle the light

# Apply the button operations
for _ in range(3):  # Three rounds of Button A
    toggle_lights(grid, button_A_positions)

# Convert the grid to a single list
final_state = [grid[i][j] for i in range(4) for j in range(4)]

# Print the final state
print(final_state)