# Initial grid
grid = [
    [1, 0],
    [1, 1]
]

# Button A controls these positions
button_a_positions = [(1, 0), (1, 1)]

# Toggle the lights controlled by Button A
for x, y in button_a_positions:
    grid[x][y] = 1 - grid[x][y]  # Toggle the light

# Flatten the grid and print the result
result = [grid[i][j] for i in range(2) for j in range(2)]
print(result)