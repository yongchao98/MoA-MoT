# Initialize a 20x20 grid with '❌'
grid = [['❌' for _ in range(20)] for _ in range(20)]

# Fill the grid according to the identified pattern
for i in range(20):
    if i == 13:
        # The row with the '❔' should have a '✅' at the 20th position
        grid[i][19] = '✅'
    else:
        # Place '✅' at the (i+6)th position
        grid[i][i+6] = '✅'

# Print the completed grid
for row in grid:
    print(''.join(row))