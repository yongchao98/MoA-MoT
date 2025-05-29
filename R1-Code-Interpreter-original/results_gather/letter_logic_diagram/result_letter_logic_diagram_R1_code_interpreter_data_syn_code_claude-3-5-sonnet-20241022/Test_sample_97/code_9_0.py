def print_solution(grid):
    for row in grid:
        print(','.join(row))

def create_initial_grid():
    return [
        ['', 'a', 'e', 'g', 'f', 'b', ''],
        ['a', '', 'g', 'f', '', '', ''],
        ['', '', '', '', 'c', 'd', 'a'],
        ['', 'f', 'b', 'c', '', 'a', ''],
        ['', 'b', 'c', '', 'a', '', 'g'],
        ['b', '', '', '', '', 'g', ''],
        ['c', '', '', 'e', 'g', '', '']
    ]

def solve():
    initial = create_initial_grid()
    grid = [[''] * 7 for _ in range(7)]
    
    # First, copy and validate initial configuration
    for i in range(7):
        for j in range(7):
            if initial[i][j] != '':
                grid[i][j] = initial[i][j]
    
    # Track available letters for each row and column
    available_row = [set('abcdefg') for _ in range(7)]
    available_col = [set('abcdefg') for _ in range(7)]
    
    # Remove used letters from available sets
    for i in range(7):
        for j in range(7):
            if grid[i][j] != '':
                if grid[i][j] in available_row[i]:
                    available_row[i].remove(grid[i][j])
                if grid[i][j] in available_col[j]:
                    available_col[j].remove(grid[i][j])
    
    def is_valid(row, col, letter):
        # Must match initial configuration
        if initial[row][col] != '' and initial[row][col] != letter:
            return False
        
        # Must be available in row and column
        if letter not in available_row[row] or letter not in available_col[col]:
            return False
        
        # Check minor diagonal constraint
        if row + col == 6:
            for i in range(7):
                j = 6 - i
                if grid[i][j] != '' and grid[i][j] != 'c':
                    return False
            if letter != 'c':
                return False
        
        return True
    
    def place_letter(row, col, letter):
        grid[row][col] = letter
        available_row[row].remove(letter)
        available_col[col].remove(letter)
    
    def remove_letter(row, col, letter):
        grid[row][col] = ''
        available_row[row].add(letter)
        available_col[col].add(letter)
    
    def find_next_cell():
        # First priority: cells on minor diagonal
        for i in range(7):
            j = 6 - i
            if grid[i][j] == '':
                return (i, j)
        
        # Second priority: cells with initial values
        for i in range(7):
            for j in range(7):
                if initial[i][j] != '' and grid[i][j] == '':
                    return (i, j)
        
        # Last priority: remaining empty cells
        for i in range(7):
            for j in range(7):
                if grid[i][j] == '':
                    return (i, j)
        return None
    
    def solve_recursive():
        cell = find_next_cell()
        if not cell:
            return True
        
        row, col = cell
        
        # If on minor diagonal, must use 'c'
        if row + col == 6:
            if is_valid(row, col, 'c'):
                place_letter(row, col, 'c')
                if solve_recursive():
                    return True
                remove_letter(row, col, 'c')
            return False
        
        # If cell has initial value, try only that value
        if initial[row][col] != '':
            letter = initial[row][col]
            if is_valid(row, col, letter):
                place_letter(row, col, letter)
                if solve_recursive():
                    return True
                remove_letter(row, col, letter)
            return False
        
        # Try each available letter
        possible_letters = available_row[row] & available_col[col]
        for letter in sorted(possible_letters):
            if is_valid(row, col, letter):
                place_letter(row, col, letter)
                if solve_recursive():
                    return True
                remove_letter(row, col, letter)
        
        return False
    
    if solve_recursive():
        print_solution(grid)
        return True
    else:
        print("No solution exists")
        return False

solve()