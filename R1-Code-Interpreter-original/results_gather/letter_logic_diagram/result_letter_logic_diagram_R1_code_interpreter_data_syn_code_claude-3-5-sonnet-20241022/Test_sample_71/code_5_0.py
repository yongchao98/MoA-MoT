def get_prefilled():
    return [
        ['', 'd', '', '', 'b', 'e', 'g'],
        ['', '', '', 'b', '', 'g', ''],
        ['f', '', '', 'e', 'g', '', ''],
        ['a', '', '', '', '', 'd', ''],
        ['', '', 'g', '', '', '', 'a'],
        ['', 'g', '', 'd', '', '', ''],
        ['', 'c', 'd', '', 'a', 'b', '']
    ]

def check_conflicts_with_prefilled(minor_letter):
    prefilled = get_prefilled()
    # Check if minor_letter conflicts with any prefilled diagonal positions
    for i in range(7):
        j = 6 - i
        if prefilled[i][j] and prefilled[i][j] != minor_letter:
            return True
    return False

def is_valid_placement(grid, row, col, letter):
    # Check row
    if letter in grid[row]:
        return False
    
    # Check column
    if letter in [grid[i][col] for i in range(7)]:
        return False
    
    # Check minor diagonal requirement
    if row + col == 6:  # If on minor diagonal
        for i in range(7):
            j = 6 - i
            if grid[i][j] and grid[i][j] != letter:
                return False
    
    return True

def solve_with_minor_letter(minor_letter):
    grid = get_prefilled()
    
    # First fill all minor diagonal positions
    for i in range(7):
        j = 6 - i
        if not grid[i][j]:
            if not is_valid_placement(grid, i, j, minor_letter):
                return None
            grid[i][j] = minor_letter
    
    def solve_remaining(pos=0):
        if pos == 49:  # 7x7 = 49 positions
            return True
            
        row, col = pos // 7, pos % 7
        
        # Skip prefilled cells and minor diagonal
        if grid[row][col] or (row + col == 6):
            return solve_remaining(pos + 1)
        
        for letter in 'abcdefg':
            if is_valid_placement(grid, row, col, letter):
                grid[row][col] = letter
                if solve_remaining(pos + 1):
                    return True
                grid[row][col] = ''
        return False
    
    if solve_remaining():
        return grid
    return None

# Try each possible letter for minor diagonal
for minor_letter in 'abcdefg':
    if not check_conflicts_with_prefilled(minor_letter):
        solution = solve_with_minor_letter(minor_letter)
        if solution:
            for row in solution:
                print(','.join(row))
            break