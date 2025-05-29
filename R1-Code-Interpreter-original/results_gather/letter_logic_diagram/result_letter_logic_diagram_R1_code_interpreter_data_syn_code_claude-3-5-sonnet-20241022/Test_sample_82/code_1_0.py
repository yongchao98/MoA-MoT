def check_prefilled_conflicts(grid):
    # Returns a set of letters that appear multiple times in any row or column
    conflicts = set()
    # Check rows
    for row in grid:
        seen = {}
        for letter in row:
            if letter != '':
                if letter in seen:
                    conflicts.add(letter)
                seen[letter] = True
    # Check columns
    for j in range(7):
        seen = {}
        for i in range(7):
            letter = grid[i][j]
            if letter != '':
                if letter in seen:
                    conflicts.add(letter)
                seen[letter] = True
    return conflicts

def is_valid_placement(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    return True

def solve(grid, diagonal_letter):
    def backtrack(pos=0):
        if pos == 49:  # 7x7 = 49 cells
            return True
            
        row, col = pos // 7, pos % 7
        
        # Skip pre-filled cells
        if grid[row][col] != '':
            return backtrack(pos + 1)
            
        # If on diagonal, must use diagonal_letter
        if row + col == 6:
            if is_valid_placement(grid, row, col, diagonal_letter):
                grid[row][col] = diagonal_letter
                if backtrack(pos + 1):
                    return True
                grid[row][col] = ''
            return False
            
        # Try each letter for non-diagonal positions
        for letter in 'abcdefg':
            if letter != diagonal_letter or row + col != 6:
                if is_valid_placement(grid, row, col, letter):
                    grid[row][col] = letter
                    if backtrack(pos + 1):
                        return True
                    grid[row][col] = ''
        return False
    
    return backtrack()

# Initial grid
initial_grid = [
    ['f','','e','','','',''],
    ['','e','','d','','','f'],
    ['','','','','g','f','a'],
    ['','d','b','','','a',''],
    ['d','b','','','a','e',''],
    ['b','g','','a','e','c','d'],
    ['','f','a','','','d','']
]

# Find existing diagonal letters
diagonal_letters = set()
for i in range(7):
    if initial_grid[i][6-i] != '':
        diagonal_letters.add(initial_grid[i][6-i])

# Try each possible diagonal letter, starting with any pre-existing ones
all_letters = list(diagonal_letters) + [l for l in 'abcdefg' if l not in diagonal_letters]

for diagonal_letter in all_letters:
    # Create a copy of the initial grid
    grid = [row[:] for row in initial_grid]
    
    # Fill all diagonal positions with the chosen letter
    valid_diagonal = True
    for i in range(7):
        if grid[i][6-i] != '' and grid[i][6-i] != diagonal_letter:
            valid_diagonal = False
            break
        if grid[i][6-i] == '':
            grid[i][6-i] = diagonal_letter
            # Check if this placement creates conflicts
            if not is_valid_placement(grid, i, 6-i, diagonal_letter):
                valid_diagonal = False
                break
    
    if not valid_diagonal:
        continue
        
    if solve(grid, diagonal_letter):
        for row in grid:
            print(','.join(row))
        break