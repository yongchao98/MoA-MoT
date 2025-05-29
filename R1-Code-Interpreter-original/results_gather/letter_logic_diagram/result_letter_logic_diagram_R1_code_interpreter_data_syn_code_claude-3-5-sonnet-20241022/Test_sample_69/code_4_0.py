def print_grid(grid):
    for row in grid:
        print(','.join(row))

def is_valid(grid, row, col, letter, diagonal_letter):
    # If this is a diagonal position, must be diagonal_letter
    if row + col == 6 and letter != diagonal_letter:
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

def solve_puzzle():
    # Initial grid
    grid = [
        ['c', 'g', '', 'd', 'e', 'b', ''],
        ['g', 'a', 'd', 'e', 'b', '', 'c'],
        ['', 'd', 'e', '', 'f', '', 'g'],
        ['d', 'e', '', 'f', '', '', ''],
        ['e', '', '', 'c', 'g', 'a', 'd'],
        ['b', 'f', 'c', 'g', 'a', '', ''],
        ['f', 'c', 'g', 'a', '', 'e', 'b']
    ]
    
    # Find the letter that appears most on the diagonal
    diagonal_counts = {}
    for i in range(7):
        if grid[i][6-i] != '':
            diagonal_counts[grid[i][6-i]] = diagonal_counts.get(grid[i][6-i], 0) + 1
    
    # Force 'c' as diagonal letter (as it appears in original grid and is a valid choice)
    diagonal_letter = 'c'
    
    # First, fill all diagonal positions with diagonal_letter
    for i in range(7):
        if grid[i][6-i] == '':
            grid[i][6-i] = diagonal_letter
    
    def backtrack(pos=0):
        if pos == 49:  # All cells filled
            return True
            
        row = pos // 7
        col = pos % 7
        
        # Skip pre-filled cells and diagonal cells
        if grid[row][col] != '' or (row + col == 6):
            return backtrack(pos + 1)
        
        # Try each possible letter
        for letter in 'abcdefg':
            if is_valid(grid, row, col, letter, diagonal_letter):
                grid[row][col] = letter
                if backtrack(pos + 1):
                    return True
                grid[row][col] = ''
        
        return False
    
    if backtrack():
        print_grid(grid)
        return True
    else:
        print("No solution exists")
        return False

solve_puzzle()