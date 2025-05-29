def print_grid(grid):
    for row in grid:
        print(','.join(cell if cell else '' for cell in row))

def check_diagonal_possible(grid, letter):
    # Check if the given letter can be used for all diagonal positions
    for i in range(7):
        if grid[i][6-i] and grid[i][6-i] != letter:
            return False
    return True

def is_valid(grid, row, col, num, diag_letter):
    # If it's a diagonal position, must match diag_letter
    if row + col == 6 and num != diag_letter:
        return False
        
    # Check row
    if num in grid[row]:
        return False
    
    # Check column
    if num in [grid[i][col] for i in range(7)]:
        return False
    
    return True

def find_empty(grid):
    for i in range(7):
        for j in range(7):
            if not grid[i][j]:
                return i, j
    return None

def solve_with_diagonal(grid, diag_letter):
    # First fill all empty diagonal positions
    temp_grid = [row[:] for row in grid]
    for i in range(7):
        if row + col == 6 and not temp_grid[i][6-i]:
            temp_grid[i][6-i] = diag_letter
    
    def solve(grid):
        empty = find_empty(grid)
        if not empty:
            return True
        
        row, col = empty
        for letter in 'abcdefg':
            if is_valid(grid, row, col, letter, diag_letter):
                grid[row][col] = letter
                if solve(grid):
                    return True
                grid[row][col] = ''
        return False
    
    if solve(temp_grid):
        return temp_grid
    return None

# Initial grid
initial_grid = [
    ['g', '', '', 'c', '', '', 'd'],
    ['', 'b', '', 'e', '', 'd', ''],
    ['b', '', 'e', 'f', 'd', '', 'a'],
    ['c', 'e', 'f', '', '', 'a', 'b'],
    ['e', 'f', 'd', 'g', 'a', 'b', ''],
    ['f', '', '', 'a', 'b', '', ''],
    ['', 'g', 'a', 'b', '', '', 'f']
]

# Try each possible letter for the diagonal
for diag_letter in 'abcdefg':
    if check_diagonal_possible(initial_grid, diag_letter):
        solution = solve_with_diagonal(initial_grid, diag_letter)
        if solution:
            print_grid(solution)
            break