def print_grid(grid):
    for row in grid:
        print(','.join(row))

def is_valid_for_diagonal(grid, letter):
    # Check if letter can be placed on all diagonal positions
    for i in range(7):
        j = 6-i
        # If position is filled, it must match
        if grid[i][j] != '' and grid[i][j] != letter:
            return False
        # Check row and column conflicts for this position
        for k in range(7):
            if k != j and grid[i][k] == letter:  # row check
                return False
            if k != i and grid[k][j] == letter:  # column check
                return False
    return True

def fill_diagonal(grid, letter):
    # Create a copy of the grid
    new_grid = [row[:] for row in grid]
    # Fill all diagonal positions
    for i in range(7):
        new_grid[i][6-i] = letter
    return new_grid

def is_valid_position(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    return True

def solve_remaining(grid):
    # Find empty position
    empty = None
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                empty = (i, j)
                break
        if empty:
            break
    
    # If no empty positions, we're done
    if not empty:
        return True
    
    row, col = empty
    # Try each letter
    for letter in 'abcdefg':
        if is_valid_position(grid, row, col, letter):
            grid[row][col] = letter
            if solve_remaining(grid):
                return True
            grid[row][col] = ''
    return False

def solve():
    # Initial grid
    initial_grid = [
        ['c', 'g', '', 'd', 'e', 'b', ''],
        ['g', 'a', 'd', 'e', 'b', '', 'c'],
        ['', 'd', 'e', '', 'f', '', 'g'],
        ['d', 'e', '', 'f', '', '', ''],
        ['e', '', '', 'c', 'g', 'a', 'd'],
        ['b', 'f', 'c', 'g', 'a', '', ''],
        ['f', 'c', 'g', 'a', '', 'e', 'b']
    ]
    
    # Try each possible letter for the diagonal
    for diagonal_letter in 'abcdefg':
        if is_valid_for_diagonal(initial_grid, diagonal_letter):
            # Create a new grid with the diagonal filled
            working_grid = fill_diagonal(initial_grid, diagonal_letter)
            
            # Try to solve the remaining positions
            if solve_remaining(working_grid):
                print_grid(working_grid)
                return True
    
    print("No solution exists")
    return False

solve()