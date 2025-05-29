def print_grid(grid):
    result = ""
    for row in grid:
        result += ",".join(row) + "\n"
    print("<<<")
    print(result.strip())
    print(">>>")

def is_valid(grid, row, col, letter, initial_grid):
    # Check if this position was pre-filled
    if initial_grid[row][col] != '' and initial_grid[row][col] != letter:
        return False
    
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Check minor diagonal (top-right to bottom-left)
    if row + col == 6:
        # Find the first letter on the minor diagonal
        diag_letter = None
        for i in range(7):
            if grid[i][6-i] != '':
                diag_letter = grid[i][6-i]
                break
        if diag_letter is None:
            # If no letter found yet, check initial grid
            for i in range(7):
                if initial_grid[i][6-i] != '':
                    diag_letter = initial_grid[i][6-i]
                    break
        if diag_letter and letter != diag_letter:
            return False
    
    return True

def solve_grid():
    # Initialize grid with initial values
    initial_grid = [
        ['a','','d','','b','g','c'],
        ['f','d','','','','','a'],
        ['d','','','','','',''],
        ['','','','','','f',''],
        ['','g','','a','f','','e'],
        ['','c','','f','','e',''],
        ['','a','','','e','','']
    ]
    
    grid = [[''] * 7 for _ in range(7)]
    
    def solve(row=0, col=0):
        if col == 7:
            row += 1
            col = 0
        if row == 7:
            return True
            
        # If cell is pre-filled, use that value
        if initial_grid[row][col] != '':
            if is_valid(grid, row, col, initial_grid[row][col], initial_grid):
                grid[row][col] = initial_grid[row][col]
                if solve(row, col + 1):
                    return True
            grid[row][col] = ''
            return False
            
        for letter in 'abcdefg':
            if is_valid(grid, row, col, letter, initial_grid):
                grid[row][col] = letter
                if solve(row, col + 1):
                    return True
                grid[row][col] = ''
        
        return False
    
    if solve():
        print_grid(grid)
        return True
    else:
        print("No solution found")
        return False

solve_grid()