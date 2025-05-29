# Initial grid setup
grid = [
    ['*', '*', '*', '*'],
    ['0', '0', '0', '*'],
    ['*', '0', '1', '*'],
    ['*', '*', '*', '*']
]

# Place the black piece at (2,1)
grid[2][1] = '0'

# Check and flip horizontally
# Row 3 (index 2)
if grid[2][1] == '0' and grid[2][3] == '0' and grid[2][2] == '1':
    grid[2][2] = '0'

# Convert grid to a single string representation
result = ','.join([cell for row in grid for cell in row])

print(result)