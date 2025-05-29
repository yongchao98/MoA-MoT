def verify_initial(grid, initial):
    return all(not initial[i][j] or grid[i][j] == initial[i][j] 
              for i in range(7) for j in range(7))

def verify_rows_cols(grid):
    # Check rows
    for row in grid:
        if sorted(row) != ['a','b','c','d','e','f','g']:
            return False
    
    # Check columns
    for j in range(7):
        col = [grid[i][j] for i in range(7)]
        if sorted(col) != ['a','b','c','d','e','f','g']:
            return False
    return True

def find_diagonal_letter(initial):
    # Find any pre-filled diagonal letters
    diag_letters = set()
    for i in range(7):
        if initial[i][6-i]:
            diag_letters.add(initial[i][6-i])
    return list(diag_letters)[0] if diag_letters else None

def solve(grid, initial, diag_letter):
    def try_fill(pos=0):
        if pos == 49:  # 7x7 = 49 cells
            return verify_initial(grid, initial) and verify_rows_cols(grid)
        
        row, col = pos // 7, pos % 7
        
        # Skip if cell is pre-filled
        if initial[row][col]:
            return try_fill(pos + 1)
        
        # If on minor diagonal, must use diag_letter
        if row + col == 6:
            if is_valid(row, col, diag_letter):
                grid[row][col] = diag_letter
                if try_fill(pos + 1):
                    return True
                grid[row][col] = ''
            return False
        
        # Try each possible letter
        for letter in 'abcdefg':
            if is_valid(row, col, letter):
                grid[row][col] = letter
                if try_fill(pos + 1):
                    return True
                grid[row][col] = ''
        return False
    
    def is_valid(row, col, letter):
        # Check row
        if letter in grid[row]:
            return False
        
        # Check column
        if letter in [grid[i][col] for i in range(7)]:
            return False
        
        return True
    
    return try_fill()

# Initial grid
initial = [
    ['b','','','c','g','f','a'],
    ['','e','c','g','f','a','b'],
    ['e','c','g','f','a','b','d'],
    ['c','','f','','b','','e'],
    ['g','f','','b','d','',''],
    ['','','b','','e','c','g'],
    ['','','d','e','','g','f']
]

# Try each possible diagonal letter
for diag_letter in 'abcdefg':
    grid = [[''] * 7 for _ in range(7)]
    # Copy initial values
    for i in range(7):
        for j in range(7):
            if initial[i][j]:
                grid[i][j] = initial[i][j]
    
    if solve(grid, initial, diag_letter):
        for row in grid:
            print(','.join(row))
        break