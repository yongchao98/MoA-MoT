def print_grid(grid):
    for row in grid:
        print(','.join(row))

def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if j != col and grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if i != row and grid[i][col] == letter:
            return False
    
    # Check minor diagonal consistency
    if row + col == 6:
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '' and grid[i][j] != letter:
                return False
    
    return True

def find_best_empty(grid):
    min_possibilities = 8
    best_pos = None
    best_letters = None
    
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                possible_letters = []
                for letter in 'abcdefg':
                    if is_valid(grid, i, j, letter):
                        possible_letters.append(letter)
                if len(possible_letters) < min_possibilities:
                    min_possibilities = len(possible_letters)
                    best_pos = (i, j)
                    best_letters = possible_letters
                if min_possibilities == 1:  # Can't get better than this
                    return best_pos, best_letters
    
    return best_pos, best_letters

def solve(grid):
    # Find cell with fewest possibilities
    pos_and_letters = find_best_empty(grid)
    if not pos_and_letters[0]:  # No empty cells left
        return True
    
    row, col = pos_and_letters[0]
    possible_letters = pos_and_letters[1]
    
    if not possible_letters:  # No valid letters for this position
        return False
    
    for letter in possible_letters:
        grid[row][col] = letter
        if solve(grid):
            return True
        grid[row][col] = ''
    
    return False

# Initialize grid
initial_grid = [
    ['', '', '', '', '', 'd', ''],
    ['', 'g', '', '', '', 'd', 'f'],
    ['', '', '', '', 'f', 'e', 'a'],
    ['b', '', '', '', 'f', 'e', ''],
    ['c', '', '', '', 'a', 'g', ''],
    ['d', 'f', '', '', 'g', '', ''],
    ['', 'e', 'a', '', 'b', '', '']
]

# Verify initial grid is valid
valid = True
for i in range(7):
    for j in range(7):
        if initial_grid[i][j] != '':
            letter = initial_grid[i][j]
            initial_grid[i][j] = ''
            if not is_valid(initial_grid, i, j, letter):
                print(f"Invalid initial position at ({i},{j}): {letter}")
                valid = False
            initial_grid[i][j] = letter

if valid:
    if solve(initial_grid):
        result = []
        for row in initial_grid:
            result.append(','.join(row))
        print('\n'.join(result))
    else:
        print("No solution exists")
else:
    print("Invalid initial grid")