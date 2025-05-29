def print_solution(grid):
    for row in grid:
        print(','.join(row))

def is_valid(grid, row, col, num):
    # Check row
    if num in grid[row]:
        return False
    
    # Check column
    if num in [grid[i][col] for i in range(7)]:
        return False
    
    # Check minor diagonal (if applicable)
    if row + col == 6:
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '' and grid[i][j] != num:
                return False
    
    return True

def matches_prefilled(grid):
    prefilled = [
        ['', 'b', 'a', '', '', '', ''],
        ['b', '', '', 'd', '', '', ''],
        ['', '', 'd', 'f', '', '', 'b'],
        ['g', 'd', '', '', 'e', '', ''],
        ['d', 'f', 'c', 'e', 'b', 'a', ''],
        ['', '', '', '', '', '', 'd'],
        ['', 'e', '', '', '', '', 'f']
    ]
    
    for i in range(7):
        for j in range(7):
            if prefilled[i][j] != '' and grid[i][j] != prefilled[i][j]:
                return False
    return True

def solve():
    # Initialize grid with prefilled values
    grid = [
        ['', 'b', 'a', '', '', '', ''],
        ['b', '', '', 'd', '', '', ''],
        ['', '', 'd', 'f', '', '', 'b'],
        ['g', 'd', '', '', 'e', '', ''],
        ['d', 'f', 'c', 'e', 'b', 'a', ''],
        ['', '', '', '', '', '', 'd'],
        ['', 'e', '', '', '', '', 'f']
    ]
    
    def backtrack(position=0):
        if position == 49:  # Reached end
            return True
        
        row = position // 7
        col = position % 7
        
        # Skip prefilled positions
        if grid[row][col] != '':
            return backtrack(position + 1)
        
        # If on minor diagonal, must match other diagonal positions
        if row + col == 6:
            diag_letter = None
            for i in range(7):
                j = 6 - i
                if grid[i][j] != '':
                    diag_letter = grid[i][j]
                    break
            
            if diag_letter:
                if is_valid(grid, row, col, diag_letter):
                    grid[row][col] = diag_letter
                    if backtrack(position + 1):
                        return True
                    grid[row][col] = ''
                return False
            else:
                # Try 'c' first for diagonal if no letter is set
                letters = ['c'] + [x for x in 'abdefg']
        else:
            letters = 'abcdefg'
        
        for letter in letters:
            if is_valid(grid, row, col, letter):
                grid[row][col] = letter
                if backtrack(position + 1):
                    return True
                grid[row][col] = ''
        
        return False
    
    # Start solving
    if backtrack():
        return grid
    return None

# Solve and print
solution = solve()
if solution and matches_prefilled(solution):
    print_solution(solution)
else:
    print("No valid solution found")