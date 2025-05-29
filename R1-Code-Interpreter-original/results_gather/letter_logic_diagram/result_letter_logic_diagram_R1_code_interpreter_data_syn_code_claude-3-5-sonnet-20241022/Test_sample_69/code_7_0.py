import numpy as np

def create_initial_grid():
    return np.array([
        ['c', 'g', '', 'd', 'e', 'b', ''],
        ['g', 'a', 'd', 'e', 'b', '', 'c'],
        ['', 'd', 'e', '', 'f', '', 'g'],
        ['d', 'e', '', 'f', '', '', ''],
        ['e', '', '', 'c', 'g', 'a', 'd'],
        ['b', 'f', 'c', 'g', 'a', '', ''],
        ['f', 'c', 'g', 'a', '', 'e', 'b']
    ])

def is_valid_diagonal(grid, letter):
    # Check if letter can be placed on all diagonal positions
    for i in range(7):
        j = 6-i
        if grid[i,j] != '' and grid[i,j] != letter:
            return False
        # Check row conflicts
        for k in range(7):
            if k != j and grid[i,k] == letter:
                return False
        # Check column conflicts
        for k in range(7):
            if k != i and grid[k,j] == letter:
                return False
    return True

def try_fill_diagonal(grid, letter):
    # Make a copy and try to fill diagonal
    new_grid = grid.copy()
    for i in range(7):
        j = 6-i
        if new_grid[i,j] == '':
            new_grid[i,j] = letter
    return new_grid

def is_valid_move(grid, row, col, num):
    # Check row
    if num in grid[row]:
        return False
    
    # Check column
    if num in grid[:,col]:
        return False
    
    return True

def find_empty(grid):
    # Find empty cell, skipping diagonal positions
    for i in range(7):
        for j in range(7):
            if grid[i,j] == '' and (i + j != 6):
                return (i, j)
    return None

def solve(grid):
    # Find empty position
    pos = find_empty(grid)
    
    # If no empty positions, we're done
    if not pos:
        return True
    
    row, col = pos
    
    # Try each letter
    for letter in 'abcdefg':
        if is_valid_move(grid, row, col, letter):
            # Make the move
            grid[row,col] = letter
            
            # Recursively try to solve rest
            if solve(grid):
                return True
                
            # If we get here, we need to backtrack
            grid[row,col] = ''
    
    return False

def solve_puzzle():
    grid = create_initial_grid()
    
    # Try each possible letter for diagonal
    for diag_letter in 'abcdefg':
        if is_valid_diagonal(grid, diag_letter):
            # Fill diagonal with this letter
            working_grid = try_fill_diagonal(grid, diag_letter)
            
            # Try to solve the rest
            if solve(working_grid):
                # Print solution
                for row in working_grid:
                    print(','.join(row))
                return True
    
    print("No solution exists")
    return False

solve_puzzle()