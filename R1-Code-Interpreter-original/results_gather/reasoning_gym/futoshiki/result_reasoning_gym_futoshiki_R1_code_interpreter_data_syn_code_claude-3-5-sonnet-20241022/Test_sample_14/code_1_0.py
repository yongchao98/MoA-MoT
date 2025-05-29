def print_grid(grid):
    for row in grid:
        print(' '.join(map(str, row)))

def is_valid(grid, row, col, num):
    # Check row
    for x in range(4):
        if grid[row][x] == num:
            return False
            
    # Check column
    for x in range(4):
        if grid[x][col] == num:
            return False
            
    # Check vertical constraint (2 > number below it)
    if row == 3 and col == 1:  # Position below 2
        if num >= 2:  # Must be less than 2
            return False
    if row == 2 and col == 1:  # Position of 2
        if grid[3][1] != 0 and num <= grid[3][1]:  # Must be greater than number below
            return False
            
    return True

def solve_futoshiki():
    # Initialize grid with known values
    grid = [
        [0, 0, 0, 1],
        [0, 0, 0, 0],
        [0, 2, 0, 4],
        [1, 0, 4, 0]
    ]
    
    def solve_recursive(row, col):
        # Move to next cell
        next_row = row + (col + 1) // 4
        next_col = (col + 1) % 4
        
        # If we've filled all cells, we're done
        if row == 4:
            return True
            
        # Skip filled cells
        if grid[row][col] != 0:
            return solve_recursive(next_row, next_col)
            
        # Try each number
        for num in range(1, 5):
            if is_valid(grid, row, col, num):
                grid[row][col] = num
                if solve_recursive(next_row, next_col):
                    return True
                grid[row][col] = 0
                
        return False
    
    # Start solving from first cell
    if solve_recursive(0, 0):
        return grid
    return None

# Solve the puzzle
solution = solve_futoshiki()

if solution:
    print("Solution found:")
    print_grid(solution)
else:
    print("No solution exists")