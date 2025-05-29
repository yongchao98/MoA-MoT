def print_grid(grid):
    for row in grid:
        print(','.join(row))

def verify_prefilled(grid, new_grid):
    for i in range(7):
        for j in range(7):
            if grid[i][j] != '' and grid[i][j] != new_grid[i][j]:
                return False
    return True

def is_valid_grid(grid):
    # Check rows
    for row in grid:
        if sorted(row) != ['a','b','c','d','e','f','g']:
            return False
    
    # Check columns
    for j in range(7):
        col = [grid[i][j] for i in range(7)]
        if sorted(col) != ['a','b','c','d','e','f','g']:
            return False
    
    # Check minor diagonal (must all be 'a')
    for i in range(7):
        if grid[i][6-i] != 'a':
            return False
            
    return True

def generate_solution():
    # Start with a base grid that has 'a' on the minor diagonal
    grid = [[''] * 7 for _ in range(7)]
    
    # Fill minor diagonal with 'a'
    for i in range(7):
        grid[i][6-i] = 'a'
    
    # Fill first row (ensuring 'a' is in the right place)
    grid[0] = ['b','e','d','g','c','f','a']
    
    # Fill remaining rows systematically while maintaining constraints
    grid[1] = ['d','f','g','c','e','a','a']
    grid[2] = ['g','d','c','f','a','b','a']
    grid[3] = ['c','g','f','a','b','e','a']
    grid[4] = ['e','c','a','b','f','d','a']
    grid[5] = ['f','a','b','e','d','g','a']
    grid[6] = ['a','b','e','d','g','c','a']
    
    return grid

# Initial grid
initial_grid = [
    ['', 'e', '', 'g', 'c', 'f', 'a'],
    ['', '', '', '', '', '', 'b'],
    ['', '', '', 'f', '', '', 'e'],
    ['', '', '', 'a', '', 'e', ''],
    ['', '', '', 'b', '', 'd', ''],
    ['f', '', 'b', 'e', '', '', 'c'],
    ['', '', '', 'd', 'g', '', 'f']
]

# Generate solution
solution = generate_solution()

# Verify solution meets all constraints
if (is_valid_grid(solution) and verify_prefilled(initial_grid, solution)):
    print_grid(solution)
else:
    print("No valid solution found")