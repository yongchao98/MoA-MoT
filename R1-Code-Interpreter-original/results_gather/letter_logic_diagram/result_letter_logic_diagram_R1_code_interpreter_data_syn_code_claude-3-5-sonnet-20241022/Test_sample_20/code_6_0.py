def validate_and_get_diagonal():
    # Initial grid
    grid = [
        ['b', '', '', 'g', '', 'c', 'e'],
        ['', 'd', 'g', 'a', '', 'e', ''],
        ['d', '', 'a', 'c', '', 'b', 'f'],
        ['g', '', '', 'e', '', '', 'd'],
        ['', '', 'e', 'b', '', '', ''],
        ['', 'e', '', 'f', 'd', 'g', 'a'],
        ['e', 'b', 'f', 'd', 'g', '', 'c']
    ]
    
    # Check diagonal positions and get the required letter
    diagonal_letter = None
    for i in range(7):
        if grid[i][6-i] != '':
            if diagonal_letter is None:
                diagonal_letter = grid[i][6-i]
            elif grid[i][6-i] != diagonal_letter:
                return None  # Inconsistent diagonal
    return grid, diagonal_letter

def is_valid(grid, row, col, letter, diagonal_letter):
    # If on diagonal, must match diagonal_letter
    if row + col == 6 and letter != diagonal_letter:
        return False
        
    # Check row
    for j in range(7):
        if j != col and grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if i != row and grid[i][col] == letter:
            return False
    
    return True

def solve(grid, diagonal_letter):
    # First, fill all diagonal positions
    for i in range(7):
        j = 6-i
        if grid[i][j] == '':
            grid[i][j] = diagonal_letter
    
    def fill_remaining(row, col):
        if col >= 7:
            row += 1
            col = 0
        if row >= 7:
            return True
            
        # Skip diagonal positions
        if row + col == 6:
            return fill_remaining(row, col + 1)
            
        if grid[row][col] != '':
            return fill_remaining(row, col + 1)
            
        used_in_row = set(grid[row])
        used_in_col = set(grid[i][col] for i in range(7))
        
        for letter in 'abcdefg':
            if letter not in used_in_row and letter not in used_in_col:
                if is_valid(grid, row, col, letter, diagonal_letter):
                    grid[row][col] = letter
                    if fill_remaining(row, col + 1):
                        return True
                    grid[row][col] = ''
        return False
    
    return fill_remaining(0, 0)

# Main solving process
result = validate_and_get_diagonal()
if result is None:
    print("Invalid initial grid - inconsistent diagonal")
else:
    grid, diagonal_letter = result
    print(f"Diagonal letter must be: {diagonal_letter}")
    
    if solve(grid, diagonal_letter):
        for row in grid:
            print(','.join(row))
    else:
        print("No solution exists")