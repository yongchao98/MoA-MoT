def verify_constraints(grid, initial_grid):
    # Verify initial values
    for i in range(7):
        for j in range(7):
            if initial_grid[i][j] != '' and grid[i][j] != initial_grid[i][j]:
                return False
    
    # Verify rows and columns
    for i in range(7):
        if len(set(grid[i])) != 7:  # Check row
            return False
        col = [grid[j][i] for j in range(7)]
        if len(set(col)) != 7:  # Check column
            return False
    
    # Verify minor diagonal
    diagonal = [grid[i][6-i] for i in range(7)]
    if len(set(diagonal)) != 1:
        return False
    
    return True

def find_minor_diagonal_letter(initial_grid):
    # Find if there's a pre-filled letter on minor diagonal
    for i in range(7):
        j = 6 - i
        if initial_grid[i][j] != '':
            return initial_grid[i][j]
    return None

def solve(grid, initial_grid, diagonal_letter):
    def try_fill(pos=0):
        if pos == 49:
            return verify_constraints(grid, initial_grid)
        
        row = pos // 7
        col = pos % 7
        
        # If pre-filled, skip
        if initial_grid[row][col] != '':
            grid[row][col] = initial_grid[row][col]
            return try_fill(pos + 1)
        
        # If on minor diagonal, must use diagonal_letter
        if row + col == 6:
            if is_valid(row, col, diagonal_letter):
                grid[row][col] = diagonal_letter
                if try_fill(pos + 1):
                    return True
                grid[row][col] = ''
            return False
        
        # Try each letter
        for letter in 'abcdefg':
            if is_valid(row, col, letter):
                grid[row][col] = letter
                if try_fill(pos + 1):
                    return True
                grid[row][col] = ''
        return False
    
    def is_valid(row, col, letter):
        # Check row
        for j in range(col):
            if grid[row][j] == letter:
                return False
        
        # Check column
        for i in range(row):
            if grid[i][col] == letter:
                return False
        
        return True
    
    return try_fill()

# Initialize grid
initial_grid = [
    ['', '', 'd', 'e', 'a', '', ''],
    ['', '', '', '', '', '', ''],
    ['', '', 'a', '', '', 'f', 'g'],
    ['', '', '', 'b', '', 'g', 'd'],
    ['', 'c', '', '', '', 'd', 'e'],
    ['c', '', 'f', '', '', '', 'a'],
    ['', 'f', 'g', '', 'e', '', 'c']
]

# Try each possible letter for the minor diagonal
for diagonal_letter in 'abcdefg':
    grid = [[''] * 7 for _ in range(7)]
    if solve(grid, initial_grid, diagonal_letter):
        for row in grid:
            print(','.join(row))
        break