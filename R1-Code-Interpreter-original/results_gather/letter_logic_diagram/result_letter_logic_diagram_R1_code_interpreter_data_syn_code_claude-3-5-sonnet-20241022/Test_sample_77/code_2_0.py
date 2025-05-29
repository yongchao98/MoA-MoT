def get_minor_diagonal_letter(grid):
    # Get the letter that should be on minor diagonal
    # or return None if no letter is fixed yet
    for i in range(7):
        if grid[i][6-i] != '':
            return grid[i][6-i]
    return None

def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # If this is on minor diagonal
    if row + col == 6:
        minor_letter = get_minor_diagonal_letter(grid)
        if minor_letter and letter != minor_letter:
            return False
    
    return True

def solve_grid(grid):
    # First, determine the minor diagonal letter if any position is filled
    minor_letter = get_minor_diagonal_letter(grid)
    
    def find_empty():
        # First fill minor diagonal if letter is known
        if minor_letter:
            for i in range(7):
                if grid[i][6-i] == '':
                    return i, 6-i
        # Then fill other positions
        for i in range(7):
            for j in range(7):
                if grid[i][j] == '':
                    return i, j
        return None
    
    def solve():
        empty = find_empty()
        if not empty:
            return True
        
        row, col = empty
        letters = 'abcdefg'
        
        # If on minor diagonal and letter not determined yet
        if row + col == 6 and not minor_letter:
            # Try each letter for the entire diagonal
            for letter in letters:
                valid = True
                # Check if this letter can be placed on all empty diagonal positions
                for i in range(7):
                    if grid[i][6-i] != '' and grid[i][6-i] != letter:
                        valid = False
                        break
                    if grid[i][6-i] == '' and not is_valid(grid, i, 6-i, letter):
                        valid = False
                        break
                if valid:
                    # Fill all empty diagonal positions with this letter
                    saved_positions = []
                    for i in range(7):
                        if grid[i][6-i] == '':
                            grid[i][6-i] = letter
                            saved_positions.append((i, 6-i))
                    if solve():
                        return True
                    # Backtrack
                    for i, j in saved_positions:
                        grid[i][j] = ''
            return False
        
        # Normal position or minor diagonal with known letter
        if row + col == 6:
            letters = minor_letter
        
        for letter in letters:
            if is_valid(grid, row, col, letter):
                grid[row][col] = letter
                if solve():
                    return True
                grid[row][col] = ''
        return False
    
    return solve()

# Initialize grid
grid = [
    ['', '', 'd', 'c', 'a', '', ''],
    ['', 'd', 'c', 'a', '', 'f', ''],
    ['d', '', '', '', '', '', 'g'],
    ['', '', '', '', '', '', 'd'],
    ['', '', '', '', '', '', ''],
    ['b', 'f', '', 'g', 'd', '', 'a'],
    ['f', '', 'g', '', '', 'a', 'b']
]

if solve_grid(grid):
    result = []
    for row in grid:
        result.append(','.join(row))
    print('<<<')
    for row in result:
        print(row)
    print('>>>')
else:
    print("No solution exists")