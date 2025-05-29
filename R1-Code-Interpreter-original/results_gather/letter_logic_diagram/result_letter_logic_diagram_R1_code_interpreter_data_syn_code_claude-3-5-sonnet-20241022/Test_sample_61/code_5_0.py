def print_solution(grid):
    for row in grid:
        print(','.join(row))

def solve_puzzle():
    # Start with empty 7x7 grid
    grid = [['' for _ in range(7)] for _ in range(7)]
    
    # Fill in the given values
    givens = [
        (0, 1, 'a'), (0, 2, 'c'), (0, 5, 'e'), (0, 6, 'd'),
        (1, 0, 'a'), (1, 2, 'f'), (1, 3, 'g'), (1, 5, 'd'), (1, 6, 'b'),
        (2, 0, 'c'), (2, 1, 'f'),
        (3, 1, 'g'), (3, 2, 'e'), (3, 3, 'd'),
        (4, 1, 'e'), (4, 2, 'd'), (4, 3, 'b'),
        (6, 0, 'd'), (6, 2, 'a')
    ]
    
    for r, c, v in givens:
        grid[r][c] = v
    
    # First, fill all minor diagonal with 'd'
    for i in range(7):
        if grid[i][6-i] != '' and grid[i][6-i] != 'd':
            return None  # Invalid initial state
        grid[i][6-i] = 'd'
    
    def is_valid(row, col, val):
        # Skip check for already filled cells
        if grid[row][col] != '':
            return False
            
        # Check row
        if val in grid[row]:
            return False
            
        # Check column
        if val in [grid[i][col] for i in range(7)]:
            return False
            
        return True
    
    def find_empty():
        for i in range(7):
            for j in range(7):
                if grid[i][j] == '':
                    return (i, j)
        return None
    
    def solve():
        empty = find_empty()
        if not empty:
            return True
            
        row, col = empty
        for val in 'abcdefg':
            if is_valid(row, col, val):
                grid[row][col] = val
                if solve():
                    return True
                grid[row][col] = ''
        return False
    
    if solve():
        return grid
    return None

# Solve and print
solution = solve_puzzle()
if solution:
    print_solution(solution)
else:
    print("No solution exists")