def is_valid_grid(grid):
    # Check rows
    for row in grid:
        for letter in 'abcdefg':
            if sum(1 for x in row if x == letter) > 1:
                return False
    
    # Check columns
    for j in range(7):
        column = [grid[i][j] for i in range(7)]
        for letter in 'abcdefg':
            if sum(1 for x in column if x == letter) > 1:
                return False
    
    # Check minor diagonal
    diagonal = [grid[i][6-i] for i in range(7)]
    first = next((x for x in diagonal if x != ''), '')
    if first and not all(x == first or x == '' for x in diagonal):
        return False
    
    return True

def matches_initial(grid, initial):
    return all(
        initial[i][j] == '' or grid[i][j] == initial[i][j]
        for i in range(7)
        for j in range(7)
    )

def solve(grid, initial, row=0, col=0):
    if col == 7:
        row += 1
        col = 0
    if row == 7:
        return is_valid_grid(grid) and matches_initial(grid, initial)

    # If pre-filled, skip to next position
    if initial[row][col] != '':
        grid[row][col] = initial[row][col]
        return solve(grid, initial, row, col + 1)

    # If on minor diagonal, try to use existing diagonal value
    if row + col == 6:
        diagonal_value = None
        for i in range(7):
            if grid[i][6-i] != '':
                diagonal_value = grid[i][6-i]
                break
        if diagonal_value:
            grid[row][col] = diagonal_value
            if solve(grid, initial, row, col + 1):
                return True
            grid[row][col] = ''
            return False

    # Try each possible letter
    for letter in 'abcdefg':
        # Quick check if letter already in row or column
        if letter in [grid[row][j] for j in range(col)] or \
           letter in [grid[i][col] for i in range(row)]:
            continue
            
        grid[row][col] = letter
        if solve(grid, initial, row, col + 1):
            return True
        grid[row][col] = ''
    
    return False

initial = [
    ['', 'b', 'd', '', 'e', 'g', 'f'],
    ['', '', 'a', 'e', 'g', '', ''],
    ['d', 'a', '', '', '', 'c', ''],
    ['a', 'e', '', '', '', 'b', ''],
    ['', '', 'f', 'c', 'b', '', 'a'],
    ['', 'f', 'c', '', 'd', '', ''],
    ['f', '', 'b', '', '', '', 'g']
]

# Create empty grid
grid = [['' for _ in range(7)] for _ in range(7)]

print('<<<')
if solve(grid, initial):
    for row in grid:
        print(','.join(row))
else:
    print("No solution exists")
print('>>>')