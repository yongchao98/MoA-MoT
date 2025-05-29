# Initialize a 20x20 grid with ❌
grid = [['❌' for _ in range(20)] for _ in range(20)]

# Fill in the check marks (✅) based on the identified pattern
for i in range(20):
    grid[19 - i][i] = '✅'

# Correct the position of the question mark (❔) to a check mark (✅)
grid[1][18] = '✅'

# Print the completed grid
for row in grid:
    print(''.join(row))