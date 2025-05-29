def find_diagonal_candidates(grid):
    possible_letters = set('abcdefg')
    # Check pre-filled cells on diagonal
    for i in range(7):
        j = 6 - i
        if grid[i][j] != '':
            possible_letters = {grid[i][j]}
            break
    
    # For each possible letter, check if it conflicts with any rows/columns
    valid_letters = set()
    for letter in possible_letters:
        valid = True
        for i in range(7):
            j = 6 - i
            # Check if conflicts with pre-filled cells
            if grid[i][j] != '' and grid[i][j] != letter:
                valid = False
                break
            # Check row conflicts
            for k in range(7):
                if k != j and grid[i][k] == letter:
                    valid = False
                    break
            if not valid:
                break
            # Check column conflicts
            for k in range(7):
                if k != i and grid[k][j] == letter:
                    valid = False
                    break
            if not valid:
                break
        if valid:
            valid_letters.add(letter)
    return valid_letters

def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    return True

def solve_with_diagonal(grid, diagonal_letter):
    # First fill the diagonal
    for i in range(7):
        j = 6 - i
        if grid[i][j] == '':
            grid[i][j] = diagonal_letter
    
    def solve_remaining():
        # Find empty cell
        empty = None
        for i in range(7):
            for j in range(7):
                if grid[i][j] == '':
                    empty = (i, j)
                    break
            if empty:
                break
        
        if not empty:
            return True
        
        row, col = empty
        for letter in 'abcdefg':
            if is_valid(grid, row, col, letter):
                grid[row][col] = letter
                if solve_remaining():
                    return True
                grid[row][col] = ''
        return False
    
    if solve_remaining():
        return True
    
    # If failed, clear the diagonal
    for i in range(7):
        j = 6 - i
        if grid[i][j] == diagonal_letter:
            grid[i][j] = ''
    return False

# Initial grid
grid = [
    ['', '', '', 'e', 'd', 'c', ''],
    ['f', '', 'e', '', '', '', ''],
    ['', '', 'd', '', 'g', '', ''],
    ['e', '', 'c', '', 'b', 'f', ''],
    ['d', '', '', 'b', '', 'a', 'e'],
    ['c', 'g', 'b', '', 'a', 'e', 'd'],
    ['', 'b', 'f', '', '', 'd', '']
]

# Save initial state
initial_grid = [[cell for cell in row] for row in grid]

# Try each possible diagonal letter
candidates = find_diagonal_candidates(grid)
solved = False
for diagonal_letter in candidates:
    # Restore initial state
    for i in range(7):
        for j in range(7):
            grid[i][j] = initial_grid[i][j]
    
    if solve_with_diagonal(grid, diagonal_letter):
        solved = True
        break

if solved:
    result = ""
    for row in grid:
        result += ','.join(row) + '\n'
    print("<<<")
    print(result.strip())
    print(">>>")
else:
    print("No solution exists")