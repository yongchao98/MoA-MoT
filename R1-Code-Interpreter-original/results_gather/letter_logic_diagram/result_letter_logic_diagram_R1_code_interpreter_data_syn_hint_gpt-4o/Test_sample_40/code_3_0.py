def is_valid(grid, row, col, letter, row_sets, col_sets):
    # Check if the letter can be placed in the given row and column
    return letter not in row_sets[row] and letter not in col_sets[col]

def solve(grid, row_sets, col_sets, row=0, col=0):
    if row == 7:
        return True
    if col == 7:
        return solve(grid, row_sets, col_sets, row + 1, 0)
    if grid[row][col] != '':
        return solve(grid, row_sets, col_sets, row, col + 1)
    
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter, row_sets, col_sets):
            grid[row][col] = letter
            row_sets[row].add(letter)
            col_sets[col].add(letter)
            if solve(grid, row_sets, col_sets, row, col + 1):
                return True
            grid[row][col] = ''
            row_sets[row].remove(letter)
            col_sets[col].remove(letter)
    
    return False

def fill_minor_diagonal(grid, letter):
    for i in range(7):
        grid[i][6-i] = letter

def find_solution(grid):
    for letter in 'abcdefg':
        # Try each letter for the minor diagonal
        temp_grid = [row[:] for row in grid]
        fill_minor_diagonal(temp_grid, letter)
        
        # Initialize row and column sets
        row_sets = [set(row) - {''} for row in temp_grid]
        col_sets = [set(temp_grid[r][c] for r in range(7)) - {''} for c in range(7)]
        
        if solve(temp_grid, row_sets, col_sets):
            return temp_grid
    return None

# Initial grid setup
grid = [
    ['e', '', '', 'g', '', 'c', 'b'],
    ['d', '', '', '', '', 'b', 'a'],
    ['', '', '', 'c', '', '', ''],
    ['', '', 'c', 'b', '', 'e', 'd'],
    ['f', '', '', '', 'e', 'd', ''],
    ['', '', '', '', '', '', ''],
    ['b', 'a', '', 'd', 'g', '', '']
]

# Find the solution
solution = find_solution(grid)

# Print the solution
if solution:
    for row in solution:
        print(','.join(row))
else:
    print("No solution found")