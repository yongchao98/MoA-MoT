def initialize_grid(initial_grid):
    # Make a copy and fill all minor diagonal positions with 'b'
    grid = [row[:] for row in initial_grid]
    for i in range(7):
        j = 6-i
        if not grid[i][j]:
            grid[i][j] = 'b'
        elif grid[i][j] != 'b':
            return None  # Invalid initial state
    return grid

def get_possible_values(grid, row, col):
    if row + col == 6:  # On minor diagonal
        return {'b'} if 'b' not in grid[row] else set()
    
    used = set()
    # Check row
    used.update(grid[row])
    # Check column
    used.update(grid[i][col] for i in range(7) if grid[i][col])
    return set('abcdefg') - used

def find_best_empty(grid):
    min_possibilities = float('inf')
    best_cell = None
    
    for i in range(7):
        for j in range(7):
            if not grid[i][j]:
                possibilities = len(get_possible_values(grid, i, j))
                if possibilities < min_possibilities:
                    min_possibilities = possibilities
                    best_cell = (i, j)
    return best_cell

def solve(grid):
    empty = find_best_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    possible_values = get_possible_values(grid, row, col)
    
    for value in sorted(possible_values):  # Try values in sorted order
        grid[row][col] = value
        if solve(grid):
            return True
        grid[row][col] = ''
    
    return False

# Initial grid
initial_grid = [
    ['a', 'g', '', '', 'f', '', 'b'],
    ['', 'd', '', '', 'e', '', 'a'],
    ['d', '', 'f', 'e', '', 'a', 'g'],
    ['c', 'f', 'e', 'b', 'a', 'g', ''],
    ['', '', 'b', '', 'g', 'd', 'c'],
    ['e', 'b', 'a', 'g', 'd', 'c', ''],
    ['b', '', 'g', 'd', 'c', 'f', 'e']
]

grid = initialize_grid(initial_grid)
if grid is None:
    print("Invalid initial configuration")
else:
    if solve(grid):
        print('<<<')
        for row in grid:
            print(','.join(row))
        print('>>>')
    else:
        print("No solution exists")