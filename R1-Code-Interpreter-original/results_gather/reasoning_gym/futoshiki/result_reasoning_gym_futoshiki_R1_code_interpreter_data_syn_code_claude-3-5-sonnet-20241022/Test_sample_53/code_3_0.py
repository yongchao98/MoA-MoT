def print_solution(grid):
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

def solve_futoshiki():
    # Initialize grid with zeros
    grid = [[0]*5 for _ in range(5)]
    
    # Place all fixed numbers first
    fixed_numbers = [
        (0, 0, 5), (0, 2, 1),
        (1, 0, 1), (1, 4, 2),
        (2, 1, 3),
        (3, 2, 3),
        (4, 1, 2), (4, 2, 5), (4, 4, 1)
    ]
    
    # Place fixed numbers
    for row, col, num in fixed_numbers:
        grid[row][col] = num
    
    def find_empty():
        for i in range(5):
            for j in range(5):
                if grid[i][j] == 0:
                    return i, j
        return None
    
    def is_fixed(row, col):
        return (row, col) in [(r, c) for r, c, _ in fixed_numbers]
    
    def backtrack():
        empty = find_empty()
        if not empty:
            return True
        
        row, col = empty
        
        # Skip if this is a fixed position
        if is_fixed(row, col):
            return backtrack()
        
        for num in range(1, 6):
            if is_valid(grid, row, col, num):
                grid[row][col] = num
                
                if backtrack():
                    return True
                
                grid[row][col] = 0
        
        return False
    
    if backtrack():
        print_solution(grid)
    else:
        print("No solution exists")

solve_futoshiki()