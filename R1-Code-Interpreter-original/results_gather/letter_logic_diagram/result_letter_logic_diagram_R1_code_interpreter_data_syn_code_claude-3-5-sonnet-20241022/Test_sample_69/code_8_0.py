def print_solution(grid):
    for row in grid:
        print(','.join(row))

def get_prefilled_diagonal(grid):
    # Return all letters that appear on the diagonal
    diagonal_letters = set()
    for i in range(7):
        if grid[i][6-i] != '':
            diagonal_letters.add(grid[i][6-i])
    return diagonal_letters

def validate_and_fill_diagonal(grid, diagonal_letter):
    # Check if diagonal_letter works with prefilled positions
    new_grid = [row[:] for row in grid]
    
    # First verify all prefilled positions
    for i in range(7):
        if new_grid[i][6-i] != '' and new_grid[i][6-i] != diagonal_letter:
            return None
    
    # Then fill all diagonal positions
    for i in range(7):
        if new_grid[i][6-i] == '':
            # Check if this placement would create conflicts
            for j in range(7):
                if j != 6-i and new_grid[i][j] == diagonal_letter:  # row check
                    return None
                if j != i and new_grid[j][6-i] == diagonal_letter:  # column check
                    return None
            new_grid[i][6-i] = diagonal_letter
            
    return new_grid

def is_valid_placement(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    return True

def find_empty(grid):
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '' and (i + j != 6):  # Skip diagonal positions
                return (i, j)
    return None

def solve_grid(grid):
    pos = find_empty(grid)
    if not pos:
        return True
    
    row, col = pos
    for letter in 'abcdefg':
        if is_valid_placement(grid, row, col, letter):
            grid[row][col] = letter
            if solve_grid(grid):
                return True
            grid[row][col] = ''
    
    return False

def solve():
    initial_grid = [
        ['c', 'g', '', 'd', 'e', 'b', ''],
        ['g', 'a', 'd', 'e', 'b', '', 'c'],
        ['', 'd', 'e', '', 'f', '', 'g'],
        ['d', 'e', '', 'f', '', '', ''],
        ['e', '', '', 'c', 'g', 'a', 'd'],
        ['b', 'f', 'c', 'g', 'a', '', ''],
        ['f', 'c', 'g', 'a', '', 'e', 'b']
    ]
    
    # Get prefilled diagonal letters
    diagonal_letters = get_prefilled_diagonal(initial_grid)
    
    # Try each possible letter for diagonal
    for diagonal_letter in 'abcdefg':
        # Skip if this letter would conflict with prefilled diagonal positions
        if diagonal_letters and diagonal_letter not in diagonal_letters:
            continue
            
        # Try to fill diagonal with this letter
        filled_grid = validate_and_fill_diagonal(initial_grid, diagonal_letter)
        if filled_grid is None:
            continue
            
        # Try to solve the rest of the grid
        if solve_grid(filled_grid):
            print_solution(filled_grid)
            return True
    
    print("No solution exists")
    return False

solve()