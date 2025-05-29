def verify_grid(grid):
    # Check rows
    for row in grid:
        if '' in row:
            continue
        if len(set(row)) != 7:
            return False
    
    # Check columns
    for j in range(7):
        col = [grid[i][j] for i in range(7)]
        if '' in col:
            continue
        if len(set(col)) != 7:
            return False
    return True

def can_use_letter_on_diagonal(grid, letter):
    # Check if the letter can be used on all diagonal positions
    # without conflicting with existing filled cells
    for i in range(7):
        j = 6 - i
        if grid[i][j] != '' and grid[i][j] != letter:
            return False
        
        # Check row and column conflicts
        for k in range(7):
            # Check row
            if k != j and grid[i][k] == letter:
                return False
            # Check column
            if k != i and grid[k][j] == letter:
                return False
    return True

def find_diagonal_letter(grid):
    initial_grid = [row[:] for row in grid]
    for letter in 'abcdefg':
        if can_use_letter_on_diagonal(initial_grid, letter):
            # Try placing this letter on all diagonal positions
            test_grid = [row[:] for row in initial_grid]
            for i in range(7):
                if test_grid[i][6-i] == '':
                    test_grid[i][6-i] = letter
            if verify_grid(test_grid):
                return letter
    return None

# Initialize grid
grid = [
    ['f', 'd', '', '', '', 'c', 'g'],
    ['d', 'e', '', 'a', 'c', '', 'f'],
    ['e', 'b', '', '', 'g', 'f', ''],
    ['', '', 'c', 'g', 'f', '', 'e'],
    ['a', 'c', 'g', 'f', 'd', 'e', 'b'],
    ['', 'g', '', 'd', 'e', '', 'a'],
    ['g', 'f', '', 'e', 'b', 'a', 'c']
]

diagonal_letter = find_diagonal_letter(grid)
if diagonal_letter:
    print(f"The letter that should appear on the diagonal is: {diagonal_letter}")
else:
    print("No valid diagonal letter found")