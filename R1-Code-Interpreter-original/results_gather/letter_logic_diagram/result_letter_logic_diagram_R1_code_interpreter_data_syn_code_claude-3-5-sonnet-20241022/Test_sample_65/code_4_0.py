def print_grid(grid):
    result = ""
    for row in grid:
        result += ",".join(row) + "\n"
    print("<<<")
    print(result.strip())
    print(">>>")

def get_minor_diagonal_letter(initial_grid):
    # Find the letter that must be on the minor diagonal
    for i in range(7):
        if initial_grid[i][6-i] != '':
            return initial_grid[i][6-i]
    return None

def is_valid(grid, row, col, letter, initial_grid, minor_diag_letter):
    # Check if this conflicts with initial values
    if initial_grid[row][col] != '' and initial_grid[row][col] != letter:
        return False

    # Check if this is on minor diagonal and must match minor_diag_letter
    if row + col == 6 and minor_diag_letter and letter != minor_diag_letter:
        return False

    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False

    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False

    return True

def solve():
    initial_grid = [
        ['a','','d','','b','g','c'],
        ['f','d','','','','','a'],
        ['d','','','','','',''],
        ['','','','','','f',''],
        ['','g','','a','f','','e'],
        ['','c','','f','','e',''],
        ['','a','','','e','','']
    ]
    
    # First, determine the minor diagonal letter
    minor_diag_letter = get_minor_diagonal_letter(initial_grid)
    if not minor_diag_letter:
        # If no letter is pre-filled on minor diagonal, we need to try each possibility
        for letter in 'abcdefg':
            if try_solve_with_minor_diag(initial_grid, letter):
                return True
    else:
        return try_solve_with_minor_diag(initial_grid, minor_diag_letter)
    
    return False

def try_solve_with_minor_diag(initial_grid, minor_diag_letter):
    grid = [[''] * 7 for _ in range(7)]
    
    # Fill in initial values and minor diagonal
    for i in range(7):
        for j in range(7):
            if initial_grid[i][j] != '':
                grid[i][j] = initial_grid[i][j]
            elif i + j == 6:  # Minor diagonal
                grid[i][j] = minor_diag_letter
    
    def backtrack(pos=0):
        if pos == 49:  # 7x7 = 49 cells
            return True
            
        row = pos // 7
        col = pos % 7
        
        # Skip pre-filled cells and minor diagonal
        if initial_grid[row][col] != '' or (row + col == 6):
            return backtrack(pos + 1)
        
        for letter in 'abcdefg':
            if is_valid(grid, row, col, letter, initial_grid, minor_diag_letter):
                grid[row][col] = letter
                if backtrack(pos + 1):
                    return True
                grid[row][col] = ''
        
        return False

    if backtrack():
        print_grid(grid)
        return True
    return False

solve()