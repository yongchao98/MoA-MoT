def print_solution(grid):
    for row in grid:
        print(','.join(row))

def is_valid_move(grid, row, col, letter):
    # Check row
    if letter in grid[row]:
        return False
    
    # Check column
    if letter in [grid[i][col] for i in range(7)]:
        return False
    
    # If on diagonal, must be 'g'
    if row + col == 6 and letter != 'g':
        return False
    
    return True

def solve_from_scratch():
    # Start with empty grid
    grid = [['' for _ in range(7)] for _ in range(7)]
    
    # Fill in the known positions from the original puzzle
    initial = [
        ['', '', 'a', 'c', '', '', 'g'],
        ['', 'a', 'c', 'b', '', 'g', ''],
        ['a', '', '', '', 'g', '', ''],
        ['c', 'b', 'f', '', 'e', 'd', 'a'],
        ['b', 'f', '', '', '', 'a', ''],
        ['', '', 'e', '', '', 'c', ''],
        ['', '', 'd', 'a', '', 'b', 'f']
    ]
    
    # First, fill all diagonal positions with 'g'
    for i in range(7):
        j = 6 - i
        grid[i][j] = 'g'
        
    # Fill in the known positions that aren't on the diagonal
    for i in range(7):
        for j in range(7):
            if initial[i][j] != '' and i + j != 6:
                grid[i][j] = initial[i][j]
    
    def solve_recursive(pos=0):
        if pos == 49:  # 7x7 = 49 positions
            return True
            
        row, col = pos // 7, pos % 7
        
        # Skip if already filled
        if grid[row][col] != '':
            return solve_recursive(pos + 1)
        
        # Try each letter
        for letter in 'abcdefg':
            if is_valid_move(grid, row, col, letter):
                grid[row][col] = letter
                if solve_recursive(pos + 1):
                    return True
                grid[row][col] = ''
                
        return False
    
    if solve_recursive():
        return grid
    return None

# Solve the puzzle
solution = solve_from_scratch()

if solution:
    print_solution(solution)
else:
    print("No solution exists")