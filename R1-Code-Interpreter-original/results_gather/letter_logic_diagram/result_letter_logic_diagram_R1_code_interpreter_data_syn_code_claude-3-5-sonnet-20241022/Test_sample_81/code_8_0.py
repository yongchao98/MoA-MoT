def is_valid_move(grid, row, col, letter, initial):
    # 1. Check if violates pre-filled position
    if initial[row][col] != '' and initial[row][col] != letter:
        return False

    # 2. Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False

    # 3. Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False

    # 4. Check minor diagonal
    if row + col == 6:
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '' and grid[i][j] != letter:
                return False

    return True

def find_next_cell(grid, initial):
    # First, fill pre-filled positions
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '' and initial[i][j] != '':
                return (i, j)
    
    # Then, fill diagonal positions
    for i in range(7):
        j = 6 - i
        if grid[i][j] == '':
            return (i, j)
    
    # Finally, fill remaining positions
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return (i, j)
    
    return None

def solve(grid, initial):
    # Find next empty cell
    pos = find_next_cell(grid, initial)
    if not pos:
        return True
    
    row, col = pos
    
    # If this is a pre-filled position, only try that letter
    if initial[row][col] != '':
        letter = initial[row][col]
        if is_valid_move(grid, row, col, letter, initial):
            grid[row][col] = letter
            if solve(grid, initial):
                return True
            grid[row][col] = ''
        return False
    
    # If on minor diagonal, try to match existing diagonal letter
    if row + col == 6:
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '':
                letter = grid[i][j]
                if is_valid_move(grid, row, col, letter, initial):
                    grid[row][col] = letter
                    if solve(grid, initial):
                        return True
                    grid[row][col] = ''
                return False
    
    # Try each possible letter
    for letter in 'abcdefg':
        if is_valid_move(grid, row, col, letter, initial):
            grid[row][col] = letter
            if solve(grid, initial):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid
initial = [
    ['', 'b', 'd', '', 'e', 'g', 'f'],
    ['', '', 'a', 'e', 'g', '', ''],
    ['d', 'a', '', '', '', 'c', ''],
    ['a', 'e', '', '', '', 'b', ''],
    ['', '', 'f', 'c', 'b', '', 'a'],
    ['', 'f', 'c', '', 'd', '', ''],
    ['f', '', 'b', '', '', '', 'g']
]

# Create working grid
grid = [['' for _ in range(7)] for _ in range(7)]

print('<<<')
if solve(grid, initial):
    for row in grid:
        print(','.join(row))
else:
    print("No solution exists")
print('>>>')