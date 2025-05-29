def print_solution(grid):
    for row in grid:
        print(','.join(row))

def get_possible_diagonal_letters(initial):
    # First, check what letters are already on the minor diagonal
    diagonal_letters = set()
    for i in range(7):
        if initial[i][6-i] != '':
            diagonal_letters.add(initial[i][6-i])
    
    if len(diagonal_letters) > 1:
        return None  # Conflict in existing diagonal
    
    # If we have a diagonal letter, that's our only option
    if len(diagonal_letters) == 1:
        return list(diagonal_letters)[0]
    
    # For each possible letter, check if it could work on diagonal
    possible_letters = []
    for letter in 'abcdefg':
        can_use = True
        for i in range(7):
            row = initial[i]
            col = [initial[r][6-i] for r in range(7) if initial[r][6-i] != '']
            # Check if letter already exists in row or column
            if letter in row or letter in col:
                if i != 6-i:  # Not on diagonal position
                    can_use = False
                    break
        if can_use:
            possible_letters.append(letter)
    
    return possible_letters

def fill_diagonal(grid, letter):
    for i in range(7):
        grid[i][6-i] = letter

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

def solve(grid):
    def find_empty():
        for i in range(7):
            for j in range(7):
                if grid[i][j] == '' and (i + j != 6):  # Skip diagonal
                    return i, j
        return None

    pos = find_empty()
    if not pos:
        return True

    row, col = pos
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid
initial = [
    ['', 'g', 'e', '', 'b', '', ''],
    ['g', 'e', 'f', 'b', '', '', ''],
    ['', '', '', '', '', '', 'g'],
    ['f', '', 'd', '', '', 'g', 'e'],
    ['', 'd', '', '', '', '', ''],
    ['d', 'c', '', '', '', '', 'b'],
    ['c', '', 'g', '', '', 'b', '']
]

# Get possible diagonal letters
diag_options = get_possible_diagonal_letters(initial)
if isinstance(diag_options, str):
    diag_options = [diag_options]

for diag_letter in diag_options:
    # Create a copy of the grid
    grid = [row[:] for row in initial]
    
    # Fill the diagonal with the chosen letter
    fill_diagonal(grid, diag_letter)
    
    # Try to solve the rest
    if solve(grid):
        print_solution(grid)
        break