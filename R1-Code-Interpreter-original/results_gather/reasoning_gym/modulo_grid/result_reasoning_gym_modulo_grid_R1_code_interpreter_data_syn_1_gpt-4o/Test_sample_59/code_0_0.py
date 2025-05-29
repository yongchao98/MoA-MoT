# Initialize a 20x20 grid with '❌'
grid = [['❌' for _ in range(20)] for _ in range(20)]

# Fill in the '✅' based on the identified pattern
for i in range(20):
    position = (i + 8) % 20  # Calculate the position of '✅' in each row
    grid[i][position] = '✅'

# Print the completed grid
for row in grid:
    print(''.join(row))