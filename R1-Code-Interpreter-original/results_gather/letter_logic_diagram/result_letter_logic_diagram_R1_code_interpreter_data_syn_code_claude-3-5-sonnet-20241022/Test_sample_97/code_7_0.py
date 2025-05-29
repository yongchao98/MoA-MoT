def print_solution(grid):
    for row in grid:
        print(','.join(row))

# Initial configuration as a constant
INITIAL = [
    ['', 'a', 'e', 'g', 'f', 'b', ''],
    ['a', '', 'g', 'f', '', '', ''],
    ['', '', '', '', 'c', 'd', 'a'],
    ['', 'f', 'b', 'c', '', 'a', ''],
    ['', 'b', 'c', '', 'a', '', 'g'],
    ['b', '', '', '', '', 'g', ''],
    ['c', '', '', 'e', 'g', '', '']
]

def solve():
    # Start with initial configuration
    grid = [['' for _ in range(7)] for _ in range(7)]
    
    # Copy initial values and track used letters in rows and columns
    rows_used = [set() for _ in range(7)]
    cols_used = [set() for _ in range(7)]
    
    # First, fill in all initial values and validate
    for i in range(7):
        for j in range(7):
            if INITIAL[i][j] != '':
                letter = INITIAL[i][j]
                if letter in rows_used[i] or letter in cols_used[j]:
                    return False
                grid[i][j] = letter
                rows_used[i].add(letter)
                cols_used[j].add(letter)
    
    def is_valid(row, col, letter):
        # Check if violates initial configuration
        if INITIAL[row][col] != '' and INITIAL[row][col] != letter:
            return False
            
        # Check row and column
        if letter in rows_used[row] or letter in cols_used[j]:
            return False
            
        # Check minor diagonal
        if row + col == 6:
            for i in range(7):
                j = 6 - i
                if grid[i][j] != '' and grid[i][j] != 'c':
                    return False
            return letter == 'c'
            
        return True
    
    def backtrack(pos=0):
        if pos == 49:
            return True
            
        row, col = pos // 7, pos % 7
        
        # Skip if cell is already filled
        if grid[row][col] != '':
            return backtrack(pos + 1)
        
        # If on minor diagonal, must use 'c'
        if row + col == 6:
            if is_valid(row, col, 'c'):
                grid[row][col] = 'c'
                rows_used[row].add('c')
                cols_used[col].add('c')
                if backtrack(pos + 1):
                    return True
                grid[row][col] = ''
                rows_used[row].remove('c')
                cols_used[col].remove('c')
            return False
        
        # Try each possible letter
        for letter in 'abcdefg':
            if is_valid(row, col, letter):
                grid[row][col] = letter
                rows_used[row].add(letter)
                cols_used[col].add(letter)
                if backtrack(pos + 1):
                    return True
                grid[row][col] = ''
                rows_used[row].remove(letter)
                cols_used[col].remove(letter)
        
        return False
    
    if backtrack():
        print_solution(grid)
        return True
    else:
        print("No solution exists")
        return False

solve()