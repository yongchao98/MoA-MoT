# Initial grid setup
grid = [
    ['*', '*', '*', '*'],
    ['*', '1', '0', '*'],
    ['0', '0', '0', '*'],
    ['*', '*', '*', '*']
]

# Place the black piece at (3,4)
grid[2][3] = '0'

# Function to print the grid in the required format
def print_grid(grid):
    result = []
    for row in grid:
        result.extend(row)
    print("<<<" + ",".join(result) + ">>>")

# Check and flip logic (no flips needed in this case)
# Print the final grid
print_grid(grid)