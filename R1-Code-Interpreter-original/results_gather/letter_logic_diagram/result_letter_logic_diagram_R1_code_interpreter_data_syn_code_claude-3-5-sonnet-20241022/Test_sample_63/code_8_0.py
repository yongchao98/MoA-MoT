def find_conflicts(grid, row, col, letter):
    conflicts = []
    # Row conflicts
    for j in range(7):
        if grid[row][j] == letter:
            conflicts.append((row, j))
    
    # Column conflicts
    for i in range(7):
        if grid[i][col] == letter:
            conflicts.append((i, col))
    
    # Diagonal conflicts
    if row + col == 6:
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '' and grid[i][j] != letter:
                conflicts.append((i, j))
    
    return conflicts

def is_complete(grid):
    return all(all(cell != '' for cell in row) for row in grid)

def solve(grid, initial_grid):
    # Step 1: Force all pre-filled positions
    for i in range(7):
        for j in range(7):
            if initial_grid[i][j] != '':
                if grid[i][j] != initial_grid[i][j]:
                    conflicts = find_conflicts(grid, i, j, initial_grid[i][j])
                    if conflicts:
                        return False
                    grid[i][j] = initial_grid[i][j]
    
    # Step 2: Find diagonal letter
    diag_letter = None
    for i in range(7):
        if grid[i][6-i] != '':
            diag_letter = grid[i][6-i]
            break
    
    if not diag_letter:
        # If no diagonal letter found, try each possible letter
        for letter in 'abcdefg':
            can_use = True
            for i in range(7):
                if grid[i][6-i] != '' and grid[i][6-i] != letter:
                    can_use = False
                    break
            if can_use:
                diag_letter = letter
                break
    
    # Step 3: Fill diagonal
    if diag_letter:
        for i in range(7):
            if grid[i][6-i] == '':
                grid[i][6-i] = diag_letter
    
    def backtrack(pos=0):
        if pos == 49:  # 7x7 = 49 positions
            return True
        
        row, col = pos // 7, pos % 7
        
        # Skip pre-filled and diagonal positions
        if grid[row][col] != '':
            return backtrack(pos + 1)
        
        # Try each letter
        for letter in 'abcdefg':
            if not find_conflicts(grid, row, col, letter):
                grid[row][col] = letter
                if backtrack(pos + 1):
                    return True
                grid[row][col] = ''
        
        return False
    
    return backtrack()

# Initialize grid
initial_grid = [
    ['', 'g', 'b', '', '', '', ''],
    ['', '', '', '', 'c', '', ''],
    ['b', '', '', '', '', 'a', ''],
    ['f', 'd', 'c', '', '', '', ''],
    ['', 'c', 'e', 'a', '', 'b', ''],
    ['', '', '', 'g', '', 'f', 'd'],
    ['e', 'a', 'g', '', 'f', '', '']
]

grid = [[''] * 7 for _ in range(7)]

if solve(grid, initial_grid):
    print("<<<")
    for row in grid:
        print(",".join(row))
    print(">>>")
else:
    print("No solution exists")