def print_grid(grid):
    for row in grid:
        print(','.join(row))

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

def check_initial_constraints(grid, diagonal_letter):
    # Check if diagonal_letter conflicts with pre-filled cells
    for i in range(7):
        j = 6-i
        if grid[i][j] != '' and grid[i][j] != diagonal_letter:
            return False
            
    # Try filling diagonal positions and check if it creates conflicts
    temp_grid = [row[:] for row in grid]
    for i in range(7):
        j = 6-i
        if temp_grid[i][j] == '':
            if not is_valid(temp_grid, i, j, diagonal_letter):
                return False
            temp_grid[i][j] = diagonal_letter
    
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
    
    def backtrack(pos=0, diagonal_letter=None):
        if pos == 49:  # All cells filled
            return True
            
        row = pos // 7
        col = pos % 7
        
        # Skip pre-filled cells
        if grid[row][col] != '':
            return backtrack(pos + 1, diagonal_letter)
        
        # If this is a diagonal position
        if row + col == 6:
            # If diagonal_letter not set yet, try each possibility
            if diagonal_letter is None:
                for letter in 'abcdefg':
                    if check_initial_constraints(grid, letter):
                        grid[row][col] = letter
                        if backtrack(pos + 1, letter):
                            return True
                        grid[row][col] = ''
            else:
                # Must use diagonal_letter
                if is_valid(grid, row, col, diagonal_letter):
                    grid[row][col] = diagonal_letter
                    if backtrack(pos + 1, diagonal_letter):
                        return True
                    grid[row][col] = ''
        else:
            # Regular position
            for letter in 'abcdefg':
                if is_valid(grid, row, col, letter):
                    grid[row][col] = letter
                    if backtrack(pos + 1, diagonal_letter):
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