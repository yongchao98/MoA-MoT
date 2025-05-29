def print_grid(grid):
    result = ""
    for row in grid:
        result += ",".join(row) + "\n"
    print("<<<")
    print(result.strip())
    print(">>>")

def is_valid(grid, row, col, letter, initial_grid):
    # First check if this conflicts with initial values
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
        # Find existing letter on minor diagonal
        diag_letter = None
        for i in range(7):
            if initial_grid[i][6-i] != '':
                diag_letter = initial_grid[i][6-i]
                break
            if grid[i][6-i] != '':
                diag_letter = grid[i][6-i]
                break
        if diag_letter and letter != diag_letter:
            return False

    return True

def copy_initial_state(grid, initial_grid):
    for i in range(7):
        for j in range(7):
            if initial_grid[i][j] != '':
                grid[i][j] = initial_grid[i][j]

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
    
    grid = [[''] * 7 for _ in range(7)]
    copy_initial_state(grid, initial_grid)
    
    def backtrack(pos=0):
        if pos == 49:  # 7x7 = 49 cells
            return True
            
        row = pos // 7
        col = pos % 7
        
        # Skip pre-filled cells
        if initial_grid[row][col] != '':
            return backtrack(pos + 1)
        
        for letter in 'abcdefg':
            if is_valid(grid, row, col, letter, initial_grid):
                grid[row][col] = letter
                if backtrack(pos + 1):
                    return True
                grid[row][col] = ''
        
        return False

    # Verify initial state is valid
    for i in range(7):
        for j in range(7):
            if initial_grid[i][j] != '':
                temp = initial_grid[i][j]
                initial_grid[i][j] = ''
                if not is_valid(initial_grid, i, j, temp, initial_grid):
                    print("Invalid initial state")
                    return False
                initial_grid[i][j] = temp

    if backtrack():
        print_grid(grid)
        return True
    else:
        print("No solution found")
        return False

solve()