def verify_prefilled(grid, initial):
    for i in range(7):
        for j in range(7):
            if initial[i][j] != '' and grid[i][j] != initial[i][j]:
                return False
    return True

def is_valid_move(grid, row, col, letter, initial):
    # Check if violates prefilled
    if initial[row][col] != '' and initial[row][col] != letter:
        return False
    
    # Check row
    for j in range(7):
        if j != col and grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if i != row and grid[i][col] == letter:
            return False
    
    # Check minor diagonal
    if row + col == 6:
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '' and grid[i][j] != letter:
                return False
    
    return True

def find_empty(grid):
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return (i, j)
    return None

def solve_puzzle():
    initial = [
        ['', 'b', 'd', '', 'e', 'g', 'f'],
        ['', '', 'a', 'e', 'g', '', ''],
        ['d', 'a', '', '', '', 'c', ''],
        ['a', 'e', '', '', '', 'b', ''],
        ['', '', 'f', 'c', 'b', '', 'a'],
        ['', 'f', 'c', '', 'd', '', ''],
        ['f', '', 'b', '', '', '', 'g']
    ]
    
    # Create working grid and copy initial values
    grid = [[''] * 7 for _ in range(7)]
    for i in range(7):
        for j in range(7):
            grid[i][j] = initial[i][j]
    
    def backtrack():
        empty = find_empty(grid)
        if not empty:
            return verify_prefilled(grid, initial)
        
        row, col = empty
        # If on minor diagonal, find the required letter
        if row + col == 6:
            for i in range(7):
                j = 6 - i
                if grid[i][j] != '':
                    letters = [grid[i][j]]
                    break
            else:
                letters = 'abcdefg'
        else:
            letters = 'abcdefg'
        
        for letter in letters:
            if is_valid_move(grid, row, col, letter, initial):
                grid[row][col] = letter
                if backtrack():
                    return True
                grid[row][col] = ''
        return False
    
    if backtrack():
        print('<<<')
        for row in grid:
            print(','.join(row))
        print('>>>')
    else:
        print('No solution exists')

solve_puzzle()