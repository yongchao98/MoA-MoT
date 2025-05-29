# Initialize the grid
grid = [
    [1, 0, 0, 0],
    [1, 1, 0, 1],
    [0, 1, 1, 1],
    [0, 1, 0, 1]
]

# Define the positions controlled by each button
button_A_positions = [(1, 2), (0, 1), (2, 1), (1, 1), (1, 3), (3, 2), (0, 2)]
button_B_positions = [(1, 2), (3, 3), (1, 3), (2, 0), (3, 1), (3, 2), (1, 1), (1, 0)]

# Function to toggle the lights
def toggle_lights(grid, positions):
    for x, y in positions:
        grid[x][y] = 1 - grid[x][y]  # Toggle the light

# Simulate the rounds
toggle_lights(grid, button_B_positions)  # Round 1: Button B
toggle_lights(grid, button_A_positions)  # Round 2: Button A
toggle_lights(grid, button_A_positions)  # Round 3: Button A

# Convert the grid to a single list
final_state = [grid[i][j] for i in range(4) for j in range(4)]

# Print the final state
print(final_state)