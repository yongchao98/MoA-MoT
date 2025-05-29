# Initialize a 20x20 grid with '❌'
grid = [['❌' for _ in range(20)] for _ in range(20)]

# Fill in the '✅' based on the identified pattern
for i in range(20):
    if i <= 10:
        grid[i][10 - i] = '✅'
    else:
        grid[i][30 - i] = '✅'

# Print the completed grid
for row in grid:
    print(''.join(row))