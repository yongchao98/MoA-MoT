def is_valid(grid, row, col, letter):
    # Check row
    if letter in [c for c in grid[row] if c]:
        return False
    
    # Check column
    if letter in [grid[i][col] for i in range(7) if grid[i][col]]:
        return False
    
    # Check minor diagonal if this cell is on it
    if row + col == 6:  # Cell is on minor diagonal
        # Find the letter that should be on minor diagonal
        diag_letter = None
        for i in range(7):
            if grid[i][6-i]:
                diag_letter = grid[i][6-i]
                break
        if diag_letter and letter != diag_letter:
            return False
    
    return True

def solve_grid(grid):
    # Find empty cell
    empty = None
    for i in range(7):
        for j in range(7):
            if not grid[i][j]:
                empty = (i, j)
                break
        if empty:
            break
    
    # If no empty cell, we're done
    if not empty:
        return True
    
    row, col = empty
    
    # Try each letter
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            
            if solve_grid(grid):
                return True
            
            grid[row][col] = ''
    
    return False

# Initialize grid
grid = [
    ['', '', 'a', '', 'c', 'f', 'b'],
    ['g', 'a', '', '', '', 'b', 'e'],
    ['', '', 'c', 'f', 'b', 'e', 'g'],
    ['d', 'c', 'f', '', 'e', 'g', 'a'],
    ['', 'f', 'b', 'e', 'g', '', ''],
    ['f', 'b', 'e', 'g', 'a', 'd', ''],
    ['b', 'e', '', 'a', 'd', '', 'f']
]

if solve_grid(grid):
    result = []
    for row in grid:
        result.append(','.join(row))
    print('\n'.join(result))
else:
    print("No solution found")