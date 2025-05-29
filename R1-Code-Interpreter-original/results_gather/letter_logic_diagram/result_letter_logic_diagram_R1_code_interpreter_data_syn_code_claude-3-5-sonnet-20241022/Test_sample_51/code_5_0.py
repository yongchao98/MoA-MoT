def verify_solution(grid, initial_grid):
    # Check if solution matches initial grid
    for i in range(7):
        for j in range(7):
            if initial_grid[i][j] != '' and grid[i][j] != initial_grid[i][j]:
                return False
    
    # Check rows
    for row in grid:
        if set(row) != set('abcdefg'):
            return False
    
    # Check columns
    for j in range(7):
        col = [grid[i][j] for i in range(7)]
        if set(col) != set('abcdefg'):
            return False
    
    # Check minor diagonal
    diagonal = [grid[i][6-i] for i in range(7)]
    if not all(x == 'd' for x in diagonal):
        return False
    
    return True

def solve(initial_grid):
    from itertools import permutations
    
    # Convert initial grid to list of lists if it's not already
    grid = [list(row) if isinstance(row, str) else row[:] for row in initial_grid]
    
    # First, fill all minor diagonal positions with 'd'
    for i in range(7):
        if grid[i][6-i] == '':
            grid[i][6-i] = 'd'
    
    # For each row, find valid permutations
    for i in range(7):
        fixed_positions = [(j, grid[i][j]) for j in range(7) if grid[i][j] != '']
        empty_positions = [j for j in range(7) if grid[i][j] == '']
        needed_letters = list(set('abcdefg') - set(x[1] for x in fixed_positions))
        
        for perm in permutations(needed_letters):
            # Try this permutation
            perm_idx = 0
            valid = True
            temp_row = grid[i][:]
            
            # Fill empty positions with permutation
            for j in empty_positions:
                temp_row[j] = perm[perm_idx]
                perm_idx += 1
            
            # Check if this creates valid columns so far
            for j in range(7):
                col = [grid[k][j] if k < i else temp_row[j] if k == i else '' for k in range(7)]
                used = [x for x in col if x != '']
                if len(used) != len(set(used)):
                    valid = False
                    break
            
            if valid:
                grid[i] = temp_row
                break
        
        if grid[i].count('') > 0:  # If we couldn't fill this row
            return None
    
    if verify_solution(grid, initial_grid):
        return grid
    return None

# Initialize grid with given values
initial_grid = [
    ['b','','','','g','e','d'],
    ['','','','','','d','b'],
    ['','','','e','','',''],
    ['','','','','b','','f'],
    ['','','d','','','f',''],
    ['e','d','','','f','',''],
    ['d','b','','f','','g','e']
]

# Solve the puzzle
solution = solve(initial_grid)

if solution:
    print('<<<')
    for row in solution:
        print(','.join(row))
    print('>>>')
else:
    print("No solution exists")