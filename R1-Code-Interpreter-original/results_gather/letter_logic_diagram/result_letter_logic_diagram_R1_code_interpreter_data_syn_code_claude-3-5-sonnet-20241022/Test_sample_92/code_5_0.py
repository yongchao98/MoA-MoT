def print_solution(grid):
    print("<<<")
    for row in grid:
        print(','.join(str(x) for x in row))
    print(">>>")

def is_valid(grid, row, col, letter, size=7):
    # Check row
    for j in range(size):
        if grid[row][j] == letter:
            return False
            
    # Check column
    for i in range(size):
        if grid[i][col] == letter:
            return False
    
    # Check minor diagonal requirement
    if row + col == size - 1:
        for i in range(size):
            if grid[i][size-1-i] and grid[i][size-1-i] != letter:
                return False
    
    return True

def find_minor_diagonal_letter(grid, size=7):
    # Find what letter is already on minor diagonal
    for i in range(size):
        if grid[i][size-1-i]:
            return grid[i][size-1-i]
    return 'a'  # If no letter found, use 'a'

def solve(grid, size=7):
    # First determine minor diagonal letter
    diagonal_letter = find_minor_diagonal_letter(grid)
    
    # Fill minor diagonal with the determined letter
    for i in range(size):
        if not grid[i][size-1-i]:
            grid[i][size-1-i] = diagonal_letter
    
    def solve_recursive(row, col):
        if col >= size:
            return solve_recursive(row + 1, 0)
        if row >= size:
            return True
            
        # Skip if cell is already filled
        if grid[row][col]:
            return solve_recursive(row, col + 1)
            
        # Try each possible letter
        for letter in 'abcdefg':
            if is_valid(grid, row, col, letter):
                grid[row][col] = letter
                if solve_recursive(row, col + 1):
                    return True
                grid[row][col] = ''
                
        return False
    
    return solve_recursive(0, 0)

# Initialize grid
initial_grid = [
    ['a', 'f', 'b', 'd', '', 'c', ''],
    ['', '', 'd', '', 'c', '', 'a'],
    ['', '', 'e', '', '', '', 'f'],
    ['', '', 'c', '', '', 'f', ''],
    ['', '', '', 'a', '', 'b', ''],
    ['', '', '', 'f', 'b', '', ''],
    ['', 'a', 'f', '', '', 'e', '']
]

# Create a copy to work with
working_grid = [row[:] for row in initial_grid]

if solve(working_grid):
    print_solution(working_grid)
else:
    # Try with different diagonal letters if first attempt fails
    for diagonal_letter in 'bcdefg':
        working_grid = [row[:] for row in initial_grid]
        # Reset diagonal cells
        for i in range(7):
            if i + (6-i) == 6:  # on minor diagonal
                if not working_grid[i][6-i]:
                    working_grid[i][6-i] = diagonal_letter
        if solve(working_grid):
            print_solution(working_grid)
            break