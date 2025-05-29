def verify_prefilled(grid, diagonal_letter):
    # Verify that the diagonal letter doesn't conflict with pre-filled cells
    for i in range(7):
        for j in range(7):
            if grid[i][j] != '':
                # Check if this pre-filled cell is on diagonal
                if i + j == 6 and grid[i][j] != diagonal_letter:
                    return False
                # Check if diagonal_letter appears in same row/column
                if grid[i][j] == diagonal_letter and i + j != 6:
                    return False
    return True

def can_place(grid, row, col, letter, diagonal_letter):
    # If on diagonal, must be diagonal_letter
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

def solve_with_diagonal(grid, diagonal_letter):
    def solve_recursive(pos=0):
        if pos == 49:
            return True
            
        row, col = pos // 7, pos % 7
        
        # Skip pre-filled cells
        if grid[row][col] != '':
            return solve_recursive(pos + 1)
        
        # If on diagonal, must use diagonal_letter
        if row + col == 6:
            if can_place(grid, row, col, diagonal_letter, diagonal_letter):
                grid[row][col] = diagonal_letter
                if solve_recursive(pos + 1):
                    return True
                grid[row][col] = ''
            return False
        
        # Try each letter for non-diagonal positions
        for letter in 'abcdefg':
            if letter != diagonal_letter and can_place(grid, row, col, letter, diagonal_letter):
                grid[row][col] = letter
                if solve_recursive(pos + 1):
                    return True
                grid[row][col] = ''
        return False
    
    return solve_recursive()

# Initial grid
grid = [
    ['e', 'a', 'c', '', '', '', ''],
    ['', '', 'g', 'f', '', 'd', 'e'],
    ['c', 'g', 'f', 'b', 'd', '', ''],
    ['', '', 'b', 'd', 'e', '', 'c'],
    ['', '', '', '', '', '', ''],
    ['', 'd', 'e', '', 'c', 'g', 'f'],
    ['d', '', '', '', 'g', 'f', '']
]

print('<<<')
solved = False
# Try each letter as the diagonal letter
for diagonal_letter in 'abcdefg':
    grid_copy = [row[:] for row in grid]
    if verify_prefilled(grid_copy, diagonal_letter):
        if solve_with_diagonal(grid_copy, diagonal_letter):
            for row in grid_copy:
                print(','.join(row))
            solved = True
            break

if not solved:
    print("No solution exists")
print('>>>')