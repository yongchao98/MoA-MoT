def verify_initial_state(grid, diagonal_letter):
    # Check if diagonal_letter conflicts with pre-filled cells
    for i in range(7):
        for j in range(7):
            if grid[i][j] != '':
                # If pre-filled cell is on diagonal, must match diagonal_letter
                if i + j == 6 and grid[i][j] != diagonal_letter:
                    return False
                # If diagonal_letter appears elsewhere, invalid
                if i + j != 6 and grid[i][j] == diagonal_letter:
                    return False
    return True

def fill_diagonal(grid, diagonal_letter):
    # Fill all diagonal positions with diagonal_letter
    for i in range(7):
        if grid[i][6-i] == '':
            grid[i][6-i] = diagonal_letter

def is_valid_move(grid, row, col, letter):
    # Skip if trying to modify diagonal
    if row + col == 6:
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

def solve_grid(grid):
    def find_empty():
        for i in range(7):
            for j in range(7):
                if grid[i][j] == '' and i + j != 6:  # Skip diagonal positions
                    return (i, j)
        return None

    pos = find_empty()
    if not pos:
        return True
        
    row, col = pos
    for letter in 'abcdefg':
        if is_valid_move(grid, row, col, letter):
            grid[row][col] = letter
            if solve_grid(grid):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid
initial_grid = [
    ['e', 'a', 'c', '', '', '', ''],
    ['', '', 'g', 'f', '', 'd', 'e'],
    ['c', 'g', 'f', 'b', 'd', '', ''],
    ['', '', 'b', 'd', 'e', '', 'c'],
    ['', '', '', '', '', '', ''],
    ['', 'd', 'e', '', 'c', 'g', 'f'],
    ['d', '', '', '', 'g', 'f', '']
]

print('<<<')
found_solution = False
for diagonal_letter in 'abcdefg':
    # Create a copy of the initial grid
    grid = [row[:] for row in initial_grid]
    
    # Verify if this diagonal_letter is possible
    if verify_initial_state(grid, diagonal_letter):
        # Fill all diagonal positions
        fill_diagonal(grid, diagonal_letter)
        
        # Try to solve the rest of the grid
        if solve_grid(grid):
            # Verify final solution
            valid = True
            for i in range(7):
                row_set = set(grid[i])
                col_set = set(grid[j][i] for j in range(7))
                if len(row_set) != 7 or len(col_set) != 7:
                    valid = False
                    break
            
            if valid:
                for row in grid:
                    print(','.join(row))
                found_solution = True
                break

if not found_solution:
    print("No solution exists")
print('>>>')