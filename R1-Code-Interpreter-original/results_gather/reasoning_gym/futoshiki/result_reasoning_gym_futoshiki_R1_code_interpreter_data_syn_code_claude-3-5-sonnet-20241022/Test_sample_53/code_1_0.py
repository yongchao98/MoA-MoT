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

def get_valid_numbers(grid, row, col):
    used = set()
    # Check row
    for j in range(5):
        if grid[row][j] != 0:
            used.add(grid[row][j])
    # Check column
    for i in range(5):
        if grid[i][col] != 0:
            used.add(grid[i][col])
    return [n for n in range(1, 6) if n not in used]

def solve_puzzle():
    # Initialize grid with fixed numbers
    grid = [
        [5, 0, 1, 0, 0],
        [1, 0, 0, 0, 2],
        [0, 3, 0, 0, 0],
        [0, 0, 3, 0, 0],
        [0, 2, 5, 0, 1]
    ]
    
    def find_most_constrained():
        min_options = 6
        best_pos = None
        for i in range(5):
            for j in range(5):
                if grid[i][j] == 0:
                    valid = get_valid_numbers(grid, i, j)
                    if len(valid) < min_options:
                        min_options = len(valid)
                        best_pos = (i, j)
        return best_pos
    
    def backtrack():
        pos = find_most_constrained()
        if not pos:
            return True
        
        row, col = pos
        valid_nums = get_valid_numbers(grid, row, col)
        
        for num in valid_nums:
            grid[row][col] = num
            if backtrack():
                return True
            grid[row][col] = 0
        
        return False
    
    if backtrack():
        print_solution(grid)
    else:
        print("No solution exists")

solve_puzzle()