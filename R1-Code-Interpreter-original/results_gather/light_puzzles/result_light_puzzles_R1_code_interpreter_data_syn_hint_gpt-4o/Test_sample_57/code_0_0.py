# Initial grid setup
grid = [
    [0, 1],
    [1, 0]
]

# Button A controls these positions (1-based index)
button_A_controls = [(1, 2), (2, 2)]

# Convert 1-based index to 0-based index for Python list
button_A_controls = [(x-1, y-1) for x, y in button_A_controls]

# Toggle the lights controlled by Button A
for x, y in button_A_controls:
    grid[x][y] = 1 - grid[x][y]  # Toggle the light

# Flatten the grid to a single list
final_state = [grid[i][j] for i in range(2) for j in range(2)]

# Print the final state
print(final_state)