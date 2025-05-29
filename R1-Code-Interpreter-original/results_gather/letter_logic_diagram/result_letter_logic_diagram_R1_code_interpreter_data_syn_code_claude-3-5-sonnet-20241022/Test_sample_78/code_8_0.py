def verify_initial_state(grid):
    # Check if any pre-filled cells on minor diagonal are not 'd'
    for i in range(7):
        j = 6 - i
        if grid[i][j] != '' and grid[i][j] != 'd':
            return False
    return True

def print_solution(grid):
    for row in grid:
        print(','.join(row))

def is_valid_placement(grid, row, col, letter):
    # Check row
    for j in range(7):
        if j != col and grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if i != row and grid[i][col] == letter:
            return False
    
    # Check minor diagonal requirement
    if row + col == 6 and letter != 'd':
        return False
        
    return True

def solve_puzzle():
    # Create initial state with known values
    solution = [
        ['', 'e', '', '', '', 'g', 'd'],  # Force 'd' at (0,6)
        ['e', 'c', 'a', 'b', 'd', '', 'f'],  # Force 'd' at (1,4)
        ['', '', 'b', 'g', 'd', 'f', 'e'],  # Force 'd' at (2,4)
        ['', 'b', 'g', 'd', '', 'e', 'c'],  # Force 'd' at (3,3)
        ['b', '', 'd', 'f', '', 'c', ''],
        ['g', 'd', '', '', 'c', 'a', 'b'],
        ['d', '', 'e', '', 'a', 'b', '']
    ]
    
    def backtrack(pos=0):
        if pos == 49:  # All cells filled
            return True
            
        row = pos // 7
        col = pos % 7
        
        # Skip pre-filled cells
        if solution[row][col] != '':
            return backtrack(pos + 1)
            
        # If on minor diagonal, must be 'd'
        if row + col == 6:
            if is_valid_placement(solution, row, col, 'd'):
                solution[row][col] = 'd'
                if backtrack(pos + 1):
                    return True
                solution[row][col] = ''
            return False
            
        # Try each possible letter
        for letter in 'abcdefg':
            if is_valid_placement(solution, row, col, letter):
                solution[row][col] = letter
                if backtrack(pos + 1):
                    return True
                solution[row][col] = ''
                
        return False
    
    if verify_initial_state(solution) and backtrack():
        return solution
    return None

# Solve the puzzle
solution = solve_puzzle()
if solution:
    print_solution(solution)
else:
    print("No solution found")