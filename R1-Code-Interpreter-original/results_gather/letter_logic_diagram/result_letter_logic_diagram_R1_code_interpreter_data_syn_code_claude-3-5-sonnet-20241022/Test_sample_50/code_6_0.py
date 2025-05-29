def check_conflicts(grid, row, col, letter):
    # Returns True if there are conflicts
    
    # Check row
    if letter in grid[row]:
        return True
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return True
    
    # Check diagonal requirement
    if row + col == 6:  # if on diagonal
        for i in range(7):
            if grid[i][6-i] and grid[i][6-i] != letter:
                return True
    elif letter in [grid[i][6-i] for i in range(7) if grid[i][6-i]]:  # if letter is used in diagonal
        if row + col != 6:  # but current position is not on diagonal
            return True
            
    return False

def find_diagonal_letter(grid):
    # Find what letter is already on diagonal
    diagonal_letters = set()
    for i in range(7):
        if grid[i][6-i]:
            diagonal_letters.add(grid[i][6-i])
    if len(diagonal_letters) > 1:
        return None
    if len(diagonal_letters) == 1:
        return list(diagonal_letters)[0]
    
    # If no letter is fixed on diagonal, find possible candidates
    candidates = set('abcdefg')
    for i in range(7):
        row, col = i, 6-i
        # Remove letters that are already in the same row or column
        for j in range(7):
            if grid[row][j]:
                candidates.discard(grid[row][j])
            if grid[j][col]:
                candidates.discard(grid[j][col])
    
    return list(candidates)[0] if candidates else None

def solve_puzzle(grid):
    def solve(pos=0):
        if pos == 49:  # 7x7 = 49 cells
            return True
            
        row, col = pos // 7, pos % 7
        
        # Skip pre-filled cells
        if grid[row][col]:
            return solve(pos + 1)
        
        # Get possible letters
        if row + col == 6:  # If on diagonal
            possible_letters = [diagonal_letter]
        else:
            possible_letters = [c for c in 'abcdefg' if c != diagonal_letter]
        
        for letter in possible_letters:
            if not check_conflicts(grid, row, col, letter):
                grid[row][col] = letter
                if solve(pos + 1):
                    return True
                grid[row][col] = ''
        
        return False
    
    # First, determine the diagonal letter
    global diagonal_letter
    diagonal_letter = find_diagonal_letter(grid)
    if not diagonal_letter:
        return False
    
    # Fill in diagonal first
    for i in range(7):
        if not grid[i][6-i]:
            grid[i][6-i] = diagonal_letter
    
    # Now solve the rest
    return solve(0)

# Initialize grid
initial = [
    ['c','','','','b','',''],
    ['','f','e','b','','g',''],
    ['','','b','','','',''],
    ['','b','a','','c','d','f'],
    ['b','a','','','','',''],
    ['','g','','','f','e',''],
    ['','c','','f','','','']
]

# Create a working copy
grid = [row[:] for row in initial]

if solve_puzzle(grid):
    for row in grid:
        print(','.join(row))
else:
    print("No solution exists")