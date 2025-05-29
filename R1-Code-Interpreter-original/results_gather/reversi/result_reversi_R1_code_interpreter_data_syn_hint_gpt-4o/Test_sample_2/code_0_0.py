# Initial grid setup
grid = [
    ['*', '0', '*', '*'],
    ['0', '0', '0', '*'],
    ['*', '0', '1', '*'],
    ['*', '*', '*', '*']
]

# Convert the grid to a single string representation
result = ''.join([''.join(row) for row in grid])
print(result)