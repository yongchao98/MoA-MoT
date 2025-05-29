def check_solution(grid):
    # Check if current solution is valid
    diagonal = [grid[i][6-i] for i in range(7)]
    if len(set(diagonal)) != 1:
        return False, "Diagonal not consistent"
    
    # Check rows and columns
    for i in range(7):
        if len(set(grid[i])) != 7:
            return False, f"Row {i} has duplicates"
        col = [grid[j][i] for j in range(7)]
        if len(set(col)) != 7:
            return False, f"Column {i} has duplicates"
    
    return True, "Valid solution"

def is_valid_move(grid, row, col, letter, diagonal_letter):
    # If on diagonal, must be diagonal_letter
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

def solve(grid, fixed_cells, diagonal_letter):
    def backtrack(pos=0):
        if pos == 49:
            valid, msg = check_solution(grid)
            return valid
        
        row, col = pos // 7, pos % 7
        
        # Skip fixed cells
        if (row, col) in fixed_cells:
            return backtrack(pos + 1)
        
        letters = [diagonal_letter] if row + col == 6 else 'abcdefg'
        for letter in letters:
            if is_valid_move(grid, row, col, letter, diagonal_letter):
                grid[row][col] = letter
                if backtrack(pos + 1):
                    return True
                grid[row][col] = ''
        return False
    
    return backtrack()

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

# Get fixed cells
fixed_cells = set()
for i in range(7):
    for j in range(7):
        if initial_grid[i][j] != '':
            fixed_cells.add((i, j))

print('<<<')
found_solution = False
# Try each letter as diagonal letter
for diagonal_letter in 'abcdefg':
    grid = [row[:] for row in initial_grid]
    
    # Check if this diagonal letter conflicts with fixed cells
    valid = True
    for i in range(7):
        if (i, 6-i) in fixed_cells and grid[i][6-i] != diagonal_letter:
            valid = False
            break
    
    if valid and solve(grid, fixed_cells, diagonal_letter):
        # Double check solution
        valid, msg = check_solution(grid)
        if valid:
            for row in grid:
                print(','.join(row))
            found_solution = True
            break

if not found_solution:
    print("No solution exists")
print('>>>')