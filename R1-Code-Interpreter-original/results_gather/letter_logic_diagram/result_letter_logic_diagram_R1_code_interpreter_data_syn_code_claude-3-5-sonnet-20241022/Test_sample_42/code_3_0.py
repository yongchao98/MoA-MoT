def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    return True

def solve_puzzle():
    # Initialize grid with the original puzzle
    grid = [
        ['', 'b', '', '', 'g', 'd', 'c'],
        ['', 'a', 'f', 'g', '', 'c', 'e'],
        ['a', 'f', '', '', 'c', 'e', 'b'],
        ['', 'g', 'd', '', 'e', '', ''],
        ['g', 'd', 'c', 'e', 'b', '', 'f'],
        ['', 'c', 'e', '', 'a', '', ''],
        ['', '', 'b', 'a', 'f', '', '']
    ]
    
    # First, fill the minor diagonal with 'c'
    for i in range(7):
        if grid[i][6-i] == '':
            grid[i][6-i] = 'c'
    
    def backtrack(row, col):
        if col == 7:
            row += 1
            col = 0
        if row == 7:
            return True
            
        # If cell is already filled, move to next cell
        if grid[row][col] != '':
            return backtrack(row, col + 1)
            
        # Try each possible letter
        for letter in 'abcdefg':
            if is_valid(grid, row, col, letter):
                grid[row][col] = letter
                if backtrack(row, col + 1):
                    return True
                grid[row][col] = ''
        return False
    
    if backtrack(0, 0):
        print('<<<')
        for row in grid:
            print(','.join(row))
        print('>>>')
        return True
    return False

solve_puzzle()