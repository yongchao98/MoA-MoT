def is_valid(grid, row, col, letter, row_sets, col_sets):
    # Check if the letter can be placed at grid[row][col]
    return letter not in row_sets[row] and letter not in col_sets[col]

def solve(grid, row_sets, col_sets, minor_diagonal_letter, row=0, col=0):
    n = len(grid)
    
    if row == n:
        return True  # Successfully filled the grid
    
    next_row, next_col = (row, col + 1) if col + 1 < n else (row + 1, 0)
    
    if grid[row][col] != '':
        return solve(grid, row_sets, col_sets, minor_diagonal_letter, next_row, next_col)
    
    for letter in 'abcdefg':
        if (col == n - row - 1 and letter != minor_diagonal_letter) or not is_valid(grid, row, col, letter, row_sets, col_sets):
            continue
        
        grid[row][col] = letter
        row_sets[row].add(letter)
        col_sets[col].add(letter)
        
        if solve(grid, row_sets, col_sets, minor_diagonal_letter, next_row, next_col):
            return True
        
        grid[row][col] = ''
        row_sets[row].remove(letter)
        col_sets[col].remove(letter)
    
    return False

def fill_grid(grid):
    n = len(grid)
    row_sets = [set(row) - {''} for row in grid]
    col_sets = [set(grid[row][col] for row in range(n)) - {''} for col in range(n)]
    
    # Determine the letter for the minor diagonal
    minor_diagonal_indices = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    possible_letters = set('abcdefg')
    
    for i, j in minor_diagonal_indices:
        if grid[i][j] != '':
            possible_letters.intersection_update(grid[i][j])
    
    minor_diagonal_letter = possible_letters.pop()
    
    # Fill the grid using backtracking
    solve(grid, row_sets, col_sets, minor_diagonal_letter)
    
    return grid

# Initial grid setup
grid = [
    ['f', '', 'e', '', '', '', ''],
    ['', 'e', '', 'd', '', '', 'f'],
    ['', '', '', '', 'g', 'f', 'a'],
    ['', 'd', 'b', '', '', 'a', ''],
    ['d', 'b', '', '', 'a', 'e', ''],
    ['b', 'g', '', 'a', 'e', 'c', 'd'],
    ['', 'f', 'a', '', '', 'd', '']
]

# Fill the grid
filled_grid = fill_grid(grid)

# Print the filled grid
for row in filled_grid:
    print(','.join(row))