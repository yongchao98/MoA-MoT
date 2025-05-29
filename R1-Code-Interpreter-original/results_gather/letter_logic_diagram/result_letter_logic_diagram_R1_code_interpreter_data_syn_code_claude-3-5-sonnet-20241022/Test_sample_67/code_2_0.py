def find_empty(grid):
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return i, j
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
    
    # For minor diagonal positions, enforce same letter
    if row + col == 6:
        diagonal_letter = None
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '':
                diagonal_letter = grid[i][j]
                break
        if diagonal_letter and letter != diagonal_letter:
            return False
    
    return True

def solve_grid(grid):
    # First fill diagonal with a consistent letter if not already filled
    diagonal_letter = None
    for i in range(7):
        j = 6 - i
        if grid[i][j] != '':
            diagonal_letter = grid[i][j]
            break
    
    if diagonal_letter is None:
        # Try each possible letter for the diagonal
        for letter in 'abcdefg':
            # Check if this letter can be used for all diagonal positions
            can_use = True
            for i in range(7):
                j = 6 - i
                if grid[i][j] != '' and grid[i][j] != letter:
                    can_use = False
                    break
                # Check if this letter already exists in the row or column
                for k in range(7):
                    if (grid[i][k] == letter and k != j) or (grid[k][j] == letter and k != i):
                        can_use = False
                        break
            if can_use:
                diagonal_letter = letter
                break
    
    if diagonal_letter is None:
        return False
    
    # Fill all empty diagonal positions
    for i in range(7):
        j = 6 - i
        if grid[i][j] == '':
            grid[i][j] = diagonal_letter
    
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve_grid(grid):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid
grid = [
    ['', 'c', 'b', '', '', 'd', ''],
    ['c', '', 'g', 'e', 'd', 'f', 'a'],
    ['b', '', '', '', '', '', 'c'],
    ['', '', '', '', 'a', '', ''],
    ['', '', 'f', '', '', '', ''],
    ['', 'f', '', 'c', '', 'g', ''],
    ['f', '', '', 'b', 'g', '', '']
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