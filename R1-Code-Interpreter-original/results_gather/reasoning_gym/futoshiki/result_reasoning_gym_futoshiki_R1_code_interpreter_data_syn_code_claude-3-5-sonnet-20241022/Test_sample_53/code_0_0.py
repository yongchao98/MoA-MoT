def is_valid(grid, row, col, num):
    # Check row
    for x in range(5):
        if grid[row][x] == num:
            return False
    
    # Check column
    for x in range(5):
        if grid[x][col] == num:
            return False
            
    return True

def solve_futoshiki(grid):
    # Given numbers
    fixed = [
        [5, 0, 1, 0, 0],
        [1, 0, 0, 0, 2],
        [0, 3, 0, 0, 0],
        [0, 0, 3, 0, 0],
        [0, 2, 5, 0, 1]
    ]
    
    def find_empty(grid):
        for i in range(5):
            for j in range(5):
                if grid[i][j] == 0:
                    return i, j
        return None
    
    def backtrack():
        find = find_empty(grid)
        if not find:
            return True
        
        row, col = find
        
        for num in range(1, 6):
            if is_valid(grid, row, col, num):
                grid[row][col] = num
                
                # Additional constraints check
                valid = True
                # Check all filled cells
                for i in range(5):
                    for j in range(5):
                        if grid[i][j] != 0:
                            # Check if it violates any fixed numbers
                            if fixed[i][j] != 0 and grid[i][j] != fixed[i][j]:
                                valid = False
                
                if valid and backtrack():
                    return True
                
                grid[row][col] = 0
        
        return False
    
    if backtrack():
        # Print in required format
        result = ""
        for i in range(5):
            row = []
            for j in range(5):
                row.append(str(grid[i][j]))
            result += "   ".join(row) + "\n"
            if i < 4:
                result += "                 \n"
        print("<<<")
        print(result.strip())
        print(">>>")
    else:
        print("No solution exists")

# Initial grid with given numbers
initial_grid = [
    [5, 0, 1, 0, 0],
    [1, 0, 0, 0, 2],
    [0, 3, 0, 0, 0],
    [0, 0, 3, 0, 0],
    [0, 2, 5, 0, 1]
]

solve_futoshiki(initial_grid)