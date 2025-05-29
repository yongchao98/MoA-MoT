def is_valid_complete(grid):
    # Check rows and columns
    for i in range(7):
        row_set = set(grid[i])
        col_set = set(grid[j][i] for j in range(7))
        if len(row_set) != 7 or len(col_set) != 7:
            return False
    
    # Check minor diagonal
    diag = [grid[i][6-i] for i in range(7)]
    if not all(x == 'g' for x in diag):
        return False
    
    return True

def solve_grid(grid):
    def is_valid_partial(row, col, letter):
        # Check row
        for j in range(7):
            if grid[row][j] == letter:
                return False
        
        # Check column
        for i in range(7):
            if grid[i][col] == letter:
                return False
        
        # Check minor diagonal
        if row + col == 6 and letter != 'g':
            return False
            
        return True
    
    def find_next_empty():
        # First fill minor diagonal
        for i in range(7):
            if grid[i][6-i] == '' and (i, 6-i) not in filled:
                return i, 6-i
        # Then other cells
        for i in range(7):
            for j in range(7):
                if grid[i][j] == '' and (i, j) not in filled:
                    return i, j
        return None
    
    filled = set()
    
    def backtrack():
        pos = find_next_empty()
        if not pos:
            return True
            
        row, col = pos
        letters = ['g'] if row + col == 6 else 'abcdefg'
        
        for letter in letters:
            if is_valid_partial(row, col, letter):
                grid[row][col] = letter
                filled.add((row, col))
                if backtrack():
                    return True
                grid[row][col] = ''
                filled.remove((row, col))
        return False
    
    # First fill all minor diagonal with 'g'
    for i in range(7):
        if grid[i][6-i] == '':
            grid[i][6-i] = 'g'
            filled.add((i, 6-i))
    
    return backtrack()

# Initialize grid
initial_grid = [
    ['b','','e','','','',''],
    ['','','','','f','',''],
    ['','','','f','g','',''],
    ['c','','f','g','b','d',''],
    ['','f','','b','d','','c'],
    ['','g','','d','','','a'],
    ['g','','','','','','f']
]

# Create working grid
grid = [row[:] for row in initial_grid]

if solve_grid(grid):
    for row in grid:
        print(','.join(row))
else:
    print("No solution exists")