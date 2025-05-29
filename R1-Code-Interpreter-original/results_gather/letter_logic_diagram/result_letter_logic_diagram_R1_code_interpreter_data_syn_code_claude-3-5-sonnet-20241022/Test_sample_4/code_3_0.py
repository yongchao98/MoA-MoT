def print_grid(grid):
    for row in grid:
        print(','.join(row))

def is_valid(grid, row, col, letter):
    # For minor diagonal positions, MUST be 'e'
    if row + col == 6:
        return letter == 'e'
    
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    return True

def solve_diagonal_first(grid):
    # First, fill all diagonal positions with 'e'
    for i in range(7):
        if grid[i][6-i] == '':
            if not is_valid(grid, i, 6-i, 'e'):
                return False
            grid[i][6-i] = 'e'
    
    # Then solve the rest
    def find_empty():
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
        for letter in 'abcdefg':
            if is_valid(grid, row, col, letter):
                grid[row][col] = letter
                if solve():
                    return True
                grid[row][col] = ''
        return False
    
    return solve()

# Initial grid
grid = [
    ['b', '', '', '', 'g', '', 'e'],
    ['', '', 'a', 'g', '', '', ''],
    ['f', 'a', '', '', '', 'b', 'd'],
    ['', 'g', '', '', '', '', ''],
    ['g', '', '', 'b', '', '', ''],
    ['', 'e', 'b', 'd', 'f', 'a', ''],
    ['', 'b', '', '', 'a', '', '']
]

# First check if initial grid is valid
initial_valid = True
for i in range(7):
    for j in range(7):
        if grid[i][j] != '':
            # Temporarily remove the letter and check if it's valid
            letter = grid[i][j]
            grid[i][j] = ''
            if not is_valid(grid, i, j, letter):
                initial_valid = False
            grid[i][j] = letter
            if not initial_valid:
                break
    if not initial_valid:
        break

if not initial_valid:
    print("No solution exists - initial configuration is invalid")
else:
    if solve_diagonal_first(grid):
        result = []
        for row in grid:
            result.append(','.join(row))
        print('\n'.join(result))
    else:
        print("No solution exists")