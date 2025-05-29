def print_grid(grid):
    for row in grid:
        print(','.join(cell if cell else '' for cell in row))

def get_minor_diagonal(grid):
    return [grid[i][6-i] for i in range(7)]

def is_valid(grid, row, col, num):
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

def can_use_diagonal_letter(grid, letter):
    # Check if the letter can be used for all diagonal positions
    # considering existing constraints
    temp_grid = [row[:] for row in grid]
    for i in range(7):
        if temp_grid[i][6-i] and temp_grid[i][6-i] != letter:
            return False
        if not temp_grid[i][6-i]:
            if not is_valid(temp_grid, i, 6-i, letter):
                return False
            temp_grid[i][6-i] = letter
    return True

def solve_with_diagonal(grid, diag_letter):
    # First fill all diagonal positions with the chosen letter
    temp_grid = [row[:] for row in grid]
    for i in range(7):
        if not temp_grid[i][6-i]:
            if not is_valid(temp_grid, i, 6-i, diag_letter):
                return None
            temp_grid[i][6-i] = diag_letter
    
    def solve_remaining(grid):
        empty = find_empty(grid)
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
    
    if solve_remaining(temp_grid):
        return temp_grid
    return None

# Initial grid
grid = [
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
    if can_use_diagonal_letter(grid, diag_letter):
        solution = solve_with_diagonal(grid, diag_letter)
        if solution:
            print_grid(solution)
            break