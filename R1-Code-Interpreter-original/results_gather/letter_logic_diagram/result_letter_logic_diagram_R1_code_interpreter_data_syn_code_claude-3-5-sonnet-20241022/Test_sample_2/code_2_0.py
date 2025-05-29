def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # For minor diagonal cells, must match the first filled diagonal cell
    if row + col == 6:
        first_diagonal = None
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '':
                first_diagonal = grid[i][j]
                break
        if first_diagonal and letter != first_diagonal:
            return False
    
    return True

def find_minor_diagonal_letter(grid):
    # Find the letter that must be on the minor diagonal
    # by checking pre-filled cells
    for i in range(7):
        j = 6 - i
        if grid[i][j] != '':
            return grid[i][j]
    return None

def solve(grid):
    # First, determine the letter that must be on the minor diagonal
    diagonal_letter = find_minor_diagonal_letter(grid)
    if not diagonal_letter:
        # If no diagonal letter is pre-filled, try each possibility
        for letter in 'abcdefg':
            # Try filling all diagonal cells with this letter
            can_use = True
            for i in range(7):
                j = 6 - i
                if grid[i][j] != '' and grid[i][j] != letter:
                    can_use = False
                    break
            if can_use:
                diagonal_letter = letter
                break
    
    # Fill all empty diagonal cells with the diagonal letter
    for i in range(7):
        j = 6 - i
        if grid[i][j] == '':
            grid[i][j] = diagonal_letter
    
    def solve_remaining(grid):
        # Find empty cell
        empty = None
        for i in range(7):
            for j in range(7):
                if grid[i][j] == '' and (i + j != 6):  # Skip diagonal cells
                    empty = (i, j)
                    break
            if empty:
                break
        
        if not empty:
            return True
        
        row, col = empty
        for letter in 'abcdefg':
            if is_valid(grid, row, col, letter):
                grid[row][col] = letter
                if solve_remaining(grid):
                    return True
                grid[row][col] = ''
        return False
    
    return solve_remaining(grid)

# Initial grid
grid = [
    ['', 'b', '', 'f', 'g', '', 'c'],
    ['b', 'e', '', '', 'a', 'c', ''],
    ['', 'f', 'g', 'a', 'c', '', 'b'],
    ['f', '', 'a', '', 'd', '', 'e'],
    ['g', '', '', 'd', '', 'e', 'f'],
    ['a', '', '', '', '', 'f', 'g'],
    ['c', 'd', '', 'e', '', 'g', 'a']
]

if solve(grid):
    result = ""
    for row in grid:
        result += ','.join(row) + '\n'
    print(result.strip())
else:
    print("No solution exists")