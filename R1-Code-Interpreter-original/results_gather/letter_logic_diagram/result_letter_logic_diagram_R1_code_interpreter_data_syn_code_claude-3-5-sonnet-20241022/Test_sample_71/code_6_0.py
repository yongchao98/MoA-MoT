def create_initial_grid():
    return [
        ['', 'd', '', '', 'b', 'e', 'g'],
        ['', '', '', 'b', '', 'g', ''],
        ['f', '', '', 'e', 'g', '', ''],
        ['a', '', '', '', '', 'd', ''],
        ['', '', 'g', '', '', '', 'a'],
        ['', 'g', '', 'd', '', '', ''],
        ['', 'c', 'd', '', 'a', 'b', '']
    ]

def check_initial_constraints(grid, minor_letter):
    # Check if minor_letter works with pre-filled cells
    initial = create_initial_grid()
    for i in range(7):
        for j in range(7):
            if initial[i][j]:
                if i + j == 6 and initial[i][j] != minor_letter:
                    return False
                if grid[i][j] != initial[i][j]:
                    return False
    return True

def is_valid(grid, row, col, letter, minor_letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Check minor diagonal requirement
    if row + col == 6:
        return letter == minor_letter
    
    return True

def find_empty(grid):
    # First fill minor diagonal
    for i in range(7):
        j = 6 - i
        if grid[i][j] == '':
            return (i, j, True)  # True indicates minor diagonal
    
    # Then fill other cells
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return (i, j, False)  # False indicates not minor diagonal
    
    return None

def solve_puzzle():
    # Try each possible letter for minor diagonal
    for minor_letter in 'abcdefg':
        grid = create_initial_grid()
        
        def solve_with_minor_letter():
            empty = find_empty(grid)
            if not empty:
                return check_initial_constraints(grid, minor_letter)
            
            row, col, is_diagonal = empty
            
            if is_diagonal:
                # For minor diagonal, only try minor_letter
                if is_valid(grid, row, col, minor_letter, minor_letter):
                    grid[row][col] = minor_letter
                    if solve_with_minor_letter():
                        return True
                    grid[row][col] = ''
            else:
                # For non-diagonal positions, try all letters
                for letter in 'abcdefg':
                    if is_valid(grid, row, col, letter, minor_letter):
                        grid[row][col] = letter
                        if solve_with_minor_letter():
                            return True
                        grid[row][col] = ''
            
            return False
        
        if solve_with_minor_letter():
            return grid
    
    return None

# Solve and print
solution = solve_puzzle()
if solution:
    for row in solution:
        print(','.join(row))
else:
    print("No solution exists")