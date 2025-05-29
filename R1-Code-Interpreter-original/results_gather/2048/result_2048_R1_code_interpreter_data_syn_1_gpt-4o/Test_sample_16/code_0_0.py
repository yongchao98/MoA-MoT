# Initial grid
grid = [
    [16, 4],
    [4, 16]
]

# Function to move left
def move_left(grid):
    for row in grid:
        if row[0] == row[1]:
            row[0] *= 2
            row[1] = 0
        elif row[0] == 0:
            row[0], row[1] = row[1], 0
    return grid

# Apply the move sequence 'll'
grid = move_left(grid)  # First 'l'
grid = move_left(grid)  # Second 'l'

# Print the final grid
print(grid)