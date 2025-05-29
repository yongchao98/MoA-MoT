def print_grid(grid):
    for row in grid:
        print(','.join(row))

def check_initial(grid):
    initial = [
        ['', 'a', 'e', 'g', 'f', 'b', ''],
        ['a', '', 'g', 'f', '', '', ''],
        ['', '', '', '', 'c', 'd', 'a'],
        ['', 'f', 'b', 'c', '', 'a', ''],
        ['', 'b', 'c', '', 'a', '', 'g'],
        ['b', '', '', '', '', 'g', ''],
        ['c', '', '', 'e', 'g', '', '']
    ]
    for i in range(7):
        for j in range(7):
            if initial[i][j] != '' and grid[i][j] != initial[i][j]:
                return False
    return True

def get_next_cell(grid):
    initial = [
        ['', 'a', 'e', 'g', 'f', 'b', ''],
        ['a', '', 'g', 'f', '', '', ''],
        ['', '', '', '', 'c', 'd', 'a'],
        ['', 'f', 'b', 'c', '', 'a', ''],
        ['', 'b', 'c', '', 'a', '', 'g'],
        ['b', '', '', '', '', 'g', ''],
        ['c', '', '', 'e', 'g', '', '']
    ]
    # First fill cells that are part of the initial configuration
    for i in range(7):
        for j in range(7):
            if initial[i][j] != '' and grid[i][j] == '':
                return (i, j)
    # Then fill remaining cells
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                return (i, j)
    return None

def is_valid(grid, row, col, letter):
    # Check initial configuration
    initial = [
        ['', 'a', 'e', 'g', 'f', 'b', ''],
        ['a', '', 'g', 'f', '', '', ''],
        ['', '', '', '', 'c', 'd', 'a'],
        ['', 'f', 'b', 'c', '', 'a', ''],
        ['', 'b', 'c', '', 'a', '', 'g'],
        ['b', '', '', '', '', 'g', ''],
        ['c', '', '', 'e', 'g', '', '']
    ]
    if initial[row][col] != '' and initial[row][col] != letter:
        return False

    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False

    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False

    # Check minor diagonal
    if row + col == 6:
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '' and grid[i][j] != letter:
                return False

    return True

def solve():
    grid = [[''] * 7 for _ in range(7)]
    
    def backtrack():
        cell = get_next_cell(grid)
        if not cell:
            return True
        
        row, col = cell
        candidates = ['c'] if row + col == 6 else 'abcdefg'
        
        for letter in candidates:
            if is_valid(grid, row, col, letter):
                grid[row][col] = letter
                if backtrack():
                    return True
                grid[row][col] = ''
        return False

    if backtrack():
        print_grid(grid)
        return True
    return False

solve()