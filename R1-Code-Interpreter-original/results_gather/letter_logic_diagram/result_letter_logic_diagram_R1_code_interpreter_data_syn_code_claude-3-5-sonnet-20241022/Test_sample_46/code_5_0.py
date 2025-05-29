def is_valid_for_position(grid, row, col, letter, diagonal_letter):
    # If this is a diagonal position, only allow diagonal_letter
    if row + col == 6:
        return letter == diagonal_letter
    
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
            
    return True

def check_initial_diagonal(grid):
    diagonal_letters = set()
    for i in range(7):
        if grid[i][6-i] != '':
            diagonal_letters.add(grid[i][6-i])
    return list(diagonal_letters)

def solve_with_diagonal(grid, diagonal_letter):
    def find_empty():
        # First fill diagonal positions
        for i in range(7):
            if grid[i][6-i] == '':
                return i, 6-i
        # Then fill other positions
        for i in range(7):
            for j in range(7):
                if grid[i][j] == '':
                    return i, j
        return None
    
    pos = find_empty()
    if not pos:
        return True
    
    row, col = pos
    
    # If this is a diagonal position, only try diagonal_letter
    if row + col == 6:
        if is_valid_for_position(grid, row, col, diagonal_letter, diagonal_letter):
            grid[row][col] = diagonal_letter
            if solve_with_diagonal(grid, diagonal_letter):
                return True
            grid[row][col] = ''
    else:
        for letter in 'abcdefg':
            if is_valid_for_position(grid, row, col, letter, diagonal_letter):
                grid[row][col] = letter
                if solve_with_diagonal(grid, diagonal_letter):
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

# Check what letters are already on diagonal
existing_diagonal = check_initial_diagonal(initial)

# Try each possible letter for diagonal
found_solution = False
for diagonal_letter in 'abcdefg':
    # Skip if conflicts with existing diagonal letters
    if existing_diagonal and diagonal_letter not in existing_diagonal:
        continue
        
    # Create a copy of the grid
    grid = [row[:] for row in initial]
    
    # Pre-fill all diagonal positions with this letter
    for i in range(7):
        if grid[i][6-i] == '':
            grid[i][6-i] = diagonal_letter
    
    # Try to solve
    if solve_with_diagonal(grid, diagonal_letter):
        for row in grid:
            print(','.join(row))
        found_solution = True
        break

if not found_solution:
    print("No solution exists")