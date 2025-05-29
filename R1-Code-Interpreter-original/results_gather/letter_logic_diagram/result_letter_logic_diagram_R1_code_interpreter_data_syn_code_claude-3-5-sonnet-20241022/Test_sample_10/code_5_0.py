def find_diagonal_letter(grid):
    # Look at pre-filled positions on minor diagonal
    diagonal_positions = [(i, 6-i) for i in range(7)]
    filled_letters = set()
    for i, j in diagonal_positions:
        if grid[i][j] != '':
            filled_letters.add(grid[i][j])
    
    # If we have a pre-filled letter on diagonal, that must be our diagonal letter
    if len(filled_letters) == 1:
        return list(filled_letters)[0]
    elif len(filled_letters) > 1:
        return None  # Conflict in diagonal
    
    # If no pre-filled diagonal letters, prefer 'f' as it appears near diagonal
    return 'f'

def create_initial_state(grid, diagonal_letter):
    # Create a new grid with the diagonal letter fixed
    new_grid = [row[:] for row in grid]
    for i in range(7):
        if new_grid[i][6-i] == '':
            new_grid[i][6-i] = diagonal_letter
    return new_grid

def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if j != col and grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if i != row and grid[i][col] == letter:
            return False
    
    # Check diagonal consistency
    if row + col == 6:
        diagonal_letter = None
        for i in range(7):
            if grid[i][6-i] != '':
                diagonal_letter = grid[i][6-i]
                break
        if diagonal_letter and letter != diagonal_letter:
            return False
    
    return True

def solve_puzzle(grid):
    def find_empty():
        # First ensure diagonal is filled
        for i in range(7):
            if grid[i][6-i] == '':
                return (i, 6-i)
        # Then fill other cells
        for i in range(7):
            for j in range(7):
                if grid[i][j] == '':
                    return (i, j)
        return None
    
    empty = find_empty()
    if not empty:
        return True
    
    row, col = empty
    candidates = ['f'] if row + col == 6 else 'abcdefg'
    
    for letter in candidates:
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve_puzzle(grid):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid
grid = [
    ['d','','','e','b','','g'],
    ['','','e','b','','','d'],
    ['','','b','','','','c'],
    ['','','a','','d','c','f'],
    ['','a','g','d','','','e'],
    ['','','','','','',''],
    ['','','','f','','','a']
]

# Find the letter that should be on the diagonal
diagonal_letter = find_diagonal_letter(grid)
if diagonal_letter:
    # Create initial state with diagonal letter fixed
    working_grid = create_initial_state(grid, diagonal_letter)
    
    # Solve the puzzle
    if solve_puzzle(working_grid):
        result = []
        for row in working_grid:
            result.append(','.join(row))
        print('\n'.join(result))
    else:
        print("No solution exists")
else:
    print("Invalid diagonal constraints")