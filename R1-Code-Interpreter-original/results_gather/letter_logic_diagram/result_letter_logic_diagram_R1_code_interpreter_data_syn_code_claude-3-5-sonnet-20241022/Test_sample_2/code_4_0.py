def verify_initial_diagonal(grid):
    # Return the letter that must be on diagonal, or None if inconsistent
    diagonal_letter = None
    for i in range(7):
        if grid[i][6-i] != '':
            if diagonal_letter is None:
                diagonal_letter = grid[i][6-i]
            elif grid[i][6-i] != diagonal_letter:
                return None
    return diagonal_letter

def can_use_diagonal_letter(grid, letter):
    # Check if this letter can be used on diagonal without conflicts
    for i in range(7):
        # Check row constraints
        if letter in grid[i] and grid[i][6-i] == '':
            return False
        # Check column constraints
        if any(grid[j][6-i] == letter for j in range(7)) and grid[i][6-i] == '':
            return False
    return True

def solve_grid():
    initial = [
        ['', 'b', '', 'f', 'g', '', 'c'],
        ['b', 'e', '', '', 'a', 'c', ''],
        ['', 'f', 'g', 'a', 'c', '', 'b'],
        ['f', '', 'a', '', 'd', '', 'e'],
        ['g', '', '', 'd', '', 'e', 'f'],
        ['a', '', '', '', '', 'f', 'g'],
        ['c', 'd', '', 'e', '', 'g', 'a']
    ]
    
    # First, verify and find diagonal letter
    diagonal_letter = verify_initial_diagonal(initial)
    if diagonal_letter is None:
        # Try each possible letter for diagonal
        for letter in 'abcdefg':
            if can_use_diagonal_letter(initial, letter):
                diagonal_letter = letter
                break
    
    if diagonal_letter is None:
        return None
    
    # Create working grid and fill diagonal
    grid = [row[:] for row in initial]
    for i in range(7):
        if grid[i][6-i] == '':
            grid[i][6-i] = diagonal_letter
    
    def available_letters(row, col):
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
        min_options = 8
        best_pos = None
        for i in range(7):
            for j in range(7):
                if grid[i][j] == '' and (i + j != 6):  # Skip diagonal
                    options = len(available_letters(i, j))
                    if options < min_options:
                        min_options = options
                        best_pos = (i, j)
        return best_pos
    
    def solve():
        pos = find_empty()
        if not pos:
            return True
        
        row, col = pos
        for letter in available_letters(row, col):
            grid[row][col] = letter
            if solve():
                return True
            grid[row][col] = ''
        return False
    
    if solve():
        result = ""
        for row in grid:
            result += ','.join(row) + "\n"
        print(result.strip())
        return grid
    return None

solve_grid()