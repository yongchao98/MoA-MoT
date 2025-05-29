def verify_initial_state(grid):
    # Check if initial diagonal positions are consistent (must all be 'e')
    diagonal_values = set()
    for i in range(7):
        if grid[i][6-i] != '':
            diagonal_values.add(grid[i][6-i])
    if len(diagonal_values) > 1:
        return False
    return True

def create_base_grid():
    # Create grid with initial values
    grid = [
        ['b', '', '', 'g', '', 'c', 'e'],
        ['', 'd', 'g', 'a', '', 'e', ''],
        ['d', '', 'a', 'c', '', 'b', 'f'],
        ['g', '', '', 'e', '', '', 'd'],
        ['', '', 'e', 'b', '', '', ''],
        ['', 'e', '', 'f', 'd', 'g', 'a'],
        ['e', 'b', 'f', 'd', 'g', '', 'c']
    ]
    
    # Force all diagonal positions to be 'e'
    for i in range(7):
        j = 6-i
        if grid[i][j] != '' and grid[i][j] != 'e':
            return None
        grid[i][j] = 'e'
    
    return grid

def is_valid(grid, row, col, letter):
    # If position is on diagonal, only 'e' is allowed
    if row + col == 6:
        return letter == 'e'
    
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    return True

def find_empty(grid):
    # First fill non-diagonal positions
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '' and (i + j != 6):
                return i, j
    return None

def solve():
    # Create and verify initial state
    grid = create_base_grid()
    if grid is None or not verify_initial_state(grid):
        return None
    
    def backtrack():
        pos = find_empty(grid)
        if not pos:
            return True
        
        row, col = pos
        for letter in 'abcdefg':
            if is_valid(grid, row, col, letter):
                grid[row][col] = letter
                if backtrack():
                    return True
                grid[row][col] = ''
        return False
    
    # Try to solve
    if backtrack():
        return grid
    return None

# Main execution
solution = solve()
if solution:
    for row in solution:
        print(','.join(row))
else:
    print("No solution exists")