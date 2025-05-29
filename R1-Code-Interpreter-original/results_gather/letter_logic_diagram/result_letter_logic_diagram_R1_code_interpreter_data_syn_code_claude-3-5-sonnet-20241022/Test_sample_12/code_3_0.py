def print_solution(grid):
    for row in grid:
        print(','.join(row))

def is_valid(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    # Minor diagonal must be 'g'
    if row + col == 6 and letter != 'g':
        return False
        
    return True

def solve(grid):
    # First, fill all minor diagonal positions with 'g'
    for i in range(7):
        if grid[i][6-i] == '':
            if not is_valid(grid, i, 6-i, 'g'):
                return False
            grid[i][6-i] = 'g'
    
    def find_empty():
        for i in range(7):
            for j in range(7):
                if grid[i][j] == '':
                    return i, j
        return None
    
    def solve_remaining():
        pos = find_empty()
        if not pos:
            return True
            
        row, col = pos
        for letter in 'abcdefg':
            if is_valid(grid, row, col, letter):
                grid[row][col] = letter
                if solve_remaining():
                    return True
                grid[row][col] = ''
        return False
    
    return solve_remaining()

# Initialize grid
grid = [
    ['b','','e','','','',''],
    ['','','','','f','',''],
    ['','','','f','g','',''],
    ['c','','f','g','b','d',''],
    ['','f','','b','d','','c'],
    ['','g','','d','','','a'],
    ['g','','','','','','f']
]

if solve(grid):
    print_solution(grid)
else:
    print("No solution exists")