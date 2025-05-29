def is_valid_move(grid, row, col, letter, diagonal_letter):
    # If on minor diagonal, must match diagonal_letter
    if row + col == 6 and letter != diagonal_letter:
        return False
    
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
            if grid[i][j] == '':
                return i, j
    return None

def solve_grid(grid, diagonal_letter):
    empty = find_empty(grid)
    if not empty:
        return True
    
    row, col = empty
    for letter in 'abcdefg':
        if is_valid_move(grid, row, col, letter, diagonal_letter):
            grid[row][col] = letter
            if solve_grid(grid, diagonal_letter):
                return True
            grid[row][col] = ''
    
    return False

def find_possible_diagonal_letters(grid):
    used_letters = set()
    # Check what letters are already fixed on the minor diagonal
    for i in range(7):
        j = 6 - i
        if grid[i][j] != '':
            used_letters.add(grid[i][j])
    
    # For each possible diagonal letter
    for diagonal_letter in 'abcdefg':
        # Create a copy of the grid
        test_grid = [row[:] for row in grid]
        # Try to fill all empty diagonal positions with this letter
        valid = True
        for i in range(7):
            j = 6 - i
            if test_grid[i][j] == '':
                test_grid[i][j] = diagonal_letter
            elif test_grid[i][j] != diagonal_letter:
                valid = False
                break
        if valid:
            # Try to solve the grid with this diagonal letter
            if solve_grid(test_grid, diagonal_letter):
                return test_grid
    return None

# Initialize grid
initial_grid = [
    ['f', 'd', '', '', '', 'c', 'g'],
    ['d', 'e', '', 'a', 'c', '', 'f'],
    ['e', 'b', '', '', 'g', 'f', ''],
    ['', '', 'c', 'g', 'f', '', 'e'],
    ['a', 'c', 'g', 'f', 'd', 'e', 'b'],
    ['', 'g', '', 'd', 'e', '', 'a'],
    ['g', 'f', '', 'e', 'b', 'a', 'c']
]

solution = find_possible_diagonal_letters(initial_grid)
if solution:
    for row in solution:
        print(','.join(row))
else:
    print("No solution exists")