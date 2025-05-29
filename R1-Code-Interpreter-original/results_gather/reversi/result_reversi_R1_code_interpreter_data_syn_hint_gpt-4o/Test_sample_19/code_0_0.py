# Initial grid setup
grid = [
    ['*', '*', '*', '*'],
    ['*', '1', '0', '*'],
    ['0', '0', '0', '0'],
    ['*', '*', '*', '*']
]

# Convert the grid to a single string representation
result = ','.join([cell for row in grid for cell in row])
print(result)