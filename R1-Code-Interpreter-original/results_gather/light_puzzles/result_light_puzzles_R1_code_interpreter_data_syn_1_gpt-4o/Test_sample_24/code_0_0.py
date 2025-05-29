# Initial grid
grid = [
    [1, 1, 1, 0, 1],
    [0, 0, 1, 0, 1],
    [1, 0, 0, 0, 0],
    [0, 0, 0, 0, 0],
    [1, 0, 0, 1, 0]
]

# Button B controls these positions (1-based index)
button_b_positions = [
    (5, 5), (3, 4), (5, 3), (5, 4), (5, 1), (3, 5), (1, 5), (2, 5)
]

# Convert to 0-based index for Python
button_b_positions = [(x-1, y-1) for x, y in button_b_positions]

# Number of rounds
rounds = 4

# Simulate the button presses
for _ in range(rounds):
    for x, y in button_b_positions:
        grid[x][y] = 1 - grid[x][y]  # Toggle the light

# Flatten the grid to a single list
final_state = [grid[i][j] for i in range(5) for j in range(5)]

# Print the final state
print(final_state)