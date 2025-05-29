def is_valid_initial(grid):
    # Check if initial grid has any contradictions
    for i in range(7):
        for j in range(7):
            if grid[i][j] != '':
                # Check row
                for k in range(7):
                    if k != j and grid[i][k] == grid[i][j]:
                        return False
                # Check column
                for k in range(7):
                    if k != i and grid[k][j] == grid[i][j]:
                        return False
    return True

def is_valid_move(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    return True

def solve():
    # Start with empty grid
    grid = [
        ['b', '', '', 'g', '', 'c', 'e'],
        ['', 'd', 'g', 'a', '', 'e', ''],
        ['d', '', 'a', 'c', '', 'b', 'f'],
        ['g', '', '', 'e', '', '', 'd'],
        ['', '', 'e', 'b', '', '', ''],
        ['', 'e', '', 'f', 'd', 'g', 'a'],
        ['e', 'b', 'f', 'd', 'g', '', 'c']
    ]
    
    if not is_valid_initial(grid):
        return None
    
    # Create a new grid with all diagonal positions set to 'e'
    new_grid = [row[:] for row in grid]
    for i in range(7):
        j = 6-i
        if new_grid[i][j] != '' and new_grid[i][j] != 'e':
            return None  # Contradiction with diagonal requirement
        new_grid[i][j] = 'e'
    
    def fill_remaining(row, col):
        if col >= 7:
            row += 1
            col = 0
        if row >= 7:
            return True
            
        # Skip if position is already filled or is on diagonal
        if new_grid[row][col] != '' or row + col == 6:
            return fill_remaining(row, col + 1)
            
        for letter in 'abcdefg':
            if letter != 'e' or (row + col == 6):  # Only use 'e' on diagonal
                if is_valid_move(new_grid, row, col, letter):
                    new_grid[row][col] = letter
                    if fill_remaining(row, col + 1):
                        return True
                    new_grid[row][col] = ''
        return False
    
    if fill_remaining(0, 0):
        return new_grid
    return None

# Solve and print result
solution = solve()
if solution:
    for row in solution:
        print(','.join(row))
else:
    print("No solution exists")