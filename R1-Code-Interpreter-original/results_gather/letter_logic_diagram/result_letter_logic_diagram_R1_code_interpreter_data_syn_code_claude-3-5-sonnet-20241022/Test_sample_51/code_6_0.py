def is_valid_column(grid, col, height):
    # Check only up to the given height
    letters = [grid[i][col] for i in range(height + 1) if grid[i][col] != '']
    return len(letters) == len(set(letters))

def is_valid_placement(grid, row, col, letter, initial_grid):
    # Check if position was pre-filled
    if initial_grid[row][col] != '' and initial_grid[row][col] != letter:
        return False
    
    # Check row
    row_letters = [x for x in grid[row] if x != '']
    if letter in row_letters:
        return False
    
    # Check column
    col_letters = [grid[i][col] for i in range(row + 1) if grid[i][col] != '']
    if letter in col_letters:
        return False
    
    # Check minor diagonal
    if row + col == 6:  # Cell is on minor diagonal
        if letter != 'd':
            return False
    
    return True

def get_candidates(grid, row, col, initial_grid):
    if initial_grid[row][col] != '':
        return [initial_grid[row][col]]
    if row + col == 6:  # Minor diagonal
        return ['d']
    return [l for l in 'abcdefg' if is_valid_placement(grid, row, col, l, initial_grid)]

def solve_grid():
    initial_grid = [
        ['b','','','','g','e','d'],
        ['','','','','','d','b'],
        ['','','','e','','',''],
        ['','','','','b','','f'],
        ['','','d','','','f',''],
        ['e','d','','','f','',''],
        ['d','b','','f','','g','e']
    ]
    
    # Create working grid
    grid = [[''] * 7 for _ in range(7)]
    
    # Copy initial values and verify minor diagonal
    for i in range(7):
        for j in range(7):
            if initial_grid[i][j] != '':
                if i + j == 6 and initial_grid[i][j] != 'd':
                    return None
                grid[i][j] = initial_grid[i][j]
    
    # Fill minor diagonal with 'd'
    for i in range(7):
        if grid[i][6-i] == '':
            grid[i][6-i] = 'd'
    
    # Stack for backtracking
    stack = []
    row, col = 0, 0
    
    while row < 7:
        if grid[row][col] == '':
            candidates = get_candidates(grid, row, col, initial_grid)
            if candidates:
                letter = candidates[0]
                grid[row][col] = letter
                stack.append((row, col, candidates[1:]))
                col = (col + 1) % 7
                if col == 0:
                    row += 1
            else:
                if not stack:
                    return None
                row, col, candidates = stack.pop()
                grid[row][col] = ''
                if candidates:
                    letter = candidates[0]
                    grid[row][col] = letter
                    stack.append((row, col, candidates[1:]))
                    col = (col + 1) % 7
                    if col == 0:
                        row += 1
        else:
            col = (col + 1) % 7
            if col == 0:
                row += 1
    
    # Verify solution
    # Check rows
    for row in grid:
        if set(row) != set('abcdefg'):
            return None
    
    # Check columns
    for j in range(7):
        if set(grid[i][j] for i in range(7)) != set('abcdefg'):
            return None
    
    # Check minor diagonal
    if not all(grid[i][6-i] == 'd' for i in range(7)):
        return None
    
    return grid

# Solve the puzzle
solution = solve_grid()

if solution:
    print('<<<')
    for row in solution:
        print(','.join(row))
    print('>>>')
else:
    print("No solution exists")