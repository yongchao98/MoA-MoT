def is_valid(grid):
    # Check rows and columns for unique numbers
    for i in range(4):
        if len(set(grid[i])) != 4 or len(set(row[i] for row in grid)) != 4:
            return False
    # Check inequalities
    if not (grid[0][0] > grid[0][1] and grid[0][2] > grid[0][3]):
        return False
    return True

def solve_futoshiki():
    # Initial grid with known numbers
    grid = [
        [0, 0, 0, 0],
        [0, 3, 0, 4],
        [0, 4, 3, 2],
        [0, 0, 0, 0]
    ]
    
    # Possible numbers to fill
    numbers = [1, 2, 3, 4]
    
    def backtrack(row, col):
        if row == 4:
            return is_valid(grid)
        if col == 4:
            return backtrack(row + 1, 0)
        if grid[row][col] != 0:
            return backtrack(row, col + 1)
        
        for num in numbers:
            grid[row][col] = num
            if is_valid(grid) and backtrack(row, col + 1):
                return True
            grid[row][col] = 0
        return False
    
    if backtrack(0, 0):
        return grid
    return None

solution = solve_futoshiki()
if solution:
    for row in solution:
        print(row)
else:
    print("No solution found")