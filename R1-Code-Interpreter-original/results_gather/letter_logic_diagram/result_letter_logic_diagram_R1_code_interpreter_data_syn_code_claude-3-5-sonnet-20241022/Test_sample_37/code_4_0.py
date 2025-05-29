def print_solution(grid):
    print('<<<')
    for row in grid:
        print(','.join(row))
    print('>>>')

def verify_initial(grid, initial):
    for i in range(7):
        for j in range(7):
            if initial[i][j] and grid[i][j] != initial[i][j]:
                return False
    return True

def get_candidates(grid, row, col, initial):
    if initial[row][col]:
        return {initial[row][col]}
        
    used_row = set(grid[row])
    used_col = set(grid[i][col] for i in range(7))
    
    # If on minor diagonal, must match existing diagonal letter if any
    if row + col == 6:
        for i in range(7):
            if grid[i][6-i]:
                return {grid[i][6-i]}
    
    return set('abcdefg') - used_row - used_col

def solve(grid, initial, filled_cells=None):
    if filled_cells is None:
        # Start with pre-filled cells
        filled_cells = []
        for i in range(7):
            for j in range(7):
                if initial[i][j]:
                    grid[i][j] = initial[i][j]
                    filled_cells.append((i, j))
    
    if len(filled_cells) == 49:
        return verify_initial(grid, initial)
    
    # Find cell with minimum candidates
    min_candidates = float('inf')
    best_pos = None
    best_candidates = None
    
    for i in range(7):
        for j in range(7):
            if (i, j) not in filled_cells:
                candidates = get_candidates(grid, i, j, initial)
                if len(candidates) < min_candidates:
                    min_candidates = len(candidates)
                    best_pos = (i, j)
                    best_candidates = candidates
                    if min_candidates == 1:  # Can't get better than this
                        break
        if min_candidates == 1:
            break
    
    if not best_candidates:
        return False
    
    i, j = best_pos
    for candidate in sorted(best_candidates):  # Sort for deterministic behavior
        grid[i][j] = candidate
        filled_cells.append((i, j))
        if solve(grid, initial, filled_cells):
            return True
        filled_cells.pop()
        grid[i][j] = ''
    
    return False

# Initialize the puzzle
initial = [
    ['c','','g','','f','e',''],
    ['a','g','','','','d',''],
    ['','','','','d','',''],
    ['','','','d','','a','g'],
    ['','e','','c','a','g','b'],
    ['e','','','a','g','b','f'],
    ['d','','a','','b','','']
]

# Create working grid
grid = [['' for _ in range(7)] for _ in range(7)]

if solve(grid, initial):
    print_solution(grid)
else:
    print("No solution exists")