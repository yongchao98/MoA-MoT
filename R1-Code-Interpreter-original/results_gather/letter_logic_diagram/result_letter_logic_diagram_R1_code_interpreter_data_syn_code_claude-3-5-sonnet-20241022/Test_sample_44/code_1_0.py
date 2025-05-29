def solve_grid():
    # Initialize grid
    grid = [
        ['b', '', '', '', '', 'a', 'f'],
        ['', '', '', '', '', '', 'b'],
        ['', 'e', 'd', 'a', '', '', ''],
        ['e', '', '', 'f', 'b', 'g', ''],
        ['', 'a', '', 'b', 'g', '', 'e'],
        ['', '', 'b', 'g', '', '', 'd'],
        ['', '', '', '', '', 'd', '']
    ]
    
    # First, verify and fill minor diagonal with 'f'
    for i in range(7):
        j = 6 - i
        if grid[i][j] != '' and grid[i][j] != 'f':
            return "No solution exists"
        grid[i][j] = 'f'
    
    def get_available_letters(row, col):
        used = set()
        # Check row
        for j in range(7):
            if grid[row][j] != '':
                used.add(grid[row][j])
        # Check column
        for i in range(7):
            if grid[i][col] != '':
                used.add(grid[i][col])
        return [c for c in 'abcdefg' if c not in used]
    
    def find_empty():
        # Find cell with fewest available options
        min_options = float('inf')
        best_pos = None
        
        for i in range(7):
            for j in range(7):
                if grid[i][j] == '' and (i + j != 6):  # Skip diagonal cells
                    available = len(get_available_letters(i, j))
                    if available < min_options:
                        min_options = available
                        best_pos = (i, j)
        return best_pos
    
    def solve():
        pos = find_empty()
        if not pos:
            return True
        
        row, col = pos
        for letter in get_available_letters(row, col):
            grid[row][col] = letter
            if solve():
                return True
            grid[row][col] = ''
        return False
    
    if solve():
        return '\n'.join(','.join(row) for row in grid)
    return "No solution exists"

print(solve_grid())