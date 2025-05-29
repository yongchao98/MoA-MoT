def print_solution(grid):
    for row in grid:
        print(','.join(row))

def is_valid_for_position(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
    
    return True

def solve_puzzle():
    # Initialize empty grid
    grid = [['' for _ in range(7)] for _ in range(7)]
    
    # Fill in the initial values
    initial_values = [
        ('g',0,0), ('c',1,2), ('b',1,5), ('g',1,6),
        ('d',2,0), ('g',2,5), ('c',3,0), ('f',3,5),
        ('e',4,0), ('b',4,2), ('g',4,3), ('f',4,4), ('c',4,6),
        ('a',5,0), ('b',5,1), ('f',5,3),
        ('f',6,2), ('c',6,4), ('a',6,6)
    ]
    
    # First, fill all initial values
    for letter, row, col in initial_values:
        grid[row][col] = letter
    
    # Find what letter should be on diagonal
    diagonal_letter = 'e'  # We'll use 'e' as specified
    
    # Fill the entire diagonal with 'e'
    for i in range(7):
        if grid[i][6-i] == '':
            grid[i][6-i] = diagonal_letter
        elif grid[i][6-i] != diagonal_letter:
            return None  # Invalid initial state
    
    def backtrack(row, col):
        if col == 7:
            row += 1
            col = 0
        if row == 7:
            return True
            
        # Skip filled cells and diagonal positions
        if grid[row][col] != '':
            return backtrack(row, col + 1)
        
        # Try each letter
        for letter in 'abcdefg':
            if is_valid_for_position(grid, row, col, letter):
                grid[row][col] = letter
                if backtrack(row, col + 1):
                    return True
                grid[row][col] = ''
        return False
    
    if backtrack(0, 0):
        return grid
    return None

# Solve and print
solution = solve_puzzle()
if solution:
    print_solution(solution)
else:
    print("No solution found")