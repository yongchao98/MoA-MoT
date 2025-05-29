def is_valid_move(grid, row, col, num, size=7):
    # Check if number already in row
    if num in grid[row]:
        return False
    
    # Check if number already in column
    if num in [grid[i][col] for i in range(size)]:
        return False
    
    # Check if on minor diagonal - must be 'a'
    if row + col == size - 1 and num != 'a':
        return False
        
    return True

def find_empty(grid, size=7):
    for i in range(size):
        for j in range(size):
            if not grid[i][j]:
                return (i, j)
    return None

def solve_puzzle(grid, size=7):
    # Find empty location
    empty = find_empty(grid)
    
    # If no empty location, puzzle is solved
    if not empty:
        return True
    
    row, col = empty
    
    # Try digits 'a' through 'g'
    for num in 'abcdefg':
        # Check if valid
        if is_valid_move(grid, row, col, num):
            # Make tentative assignment
            grid[row][col] = num
            
            # Return if success
            if solve_puzzle(grid):
                return True
            
            # Failure, unmake & try again
            grid[row][col] = ''
            
    # Trigger backtracking
    return False

# Initialize grid with given values
initial_grid = [
    ['a', 'f', 'b', 'd', '', 'c', 'a'],
    ['', '', 'd', '', 'c', '', 'a'],
    ['', '', 'e', '', '', '', 'a'],
    ['', '', 'c', '', '', 'f', 'a'],
    ['', '', '', 'a', '', 'b', 'a'],
    ['', '', '', 'f', 'b', '', 'a'],
    ['', 'a', 'f', '', '', 'e', 'a']
]

def verify_solution(grid, size=7):
    # Check rows
    for row in grid:
        if sorted(row) != sorted('abcdefg'):
            return False
    
    # Check columns
    for j in range(size):
        col = [grid[i][j] for i in range(size)]
        if sorted(col) != sorted('abcdefg'):
            return False
    
    # Check minor diagonal
    diagonal = [grid[i][size-1-i] for i in range(size)]
    if not all(x == 'a' for x in diagonal):
        return False
    
    return True

# Solve and print
if solve_puzzle(initial_grid):
    print("<<<")
    for row in initial_grid:
        print(','.join(row))
    print(">>>")
    
    # Verify solution
    if not verify_solution(initial_grid):
        print("Solution verification failed!")
else:
    print("No solution exists")