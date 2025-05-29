def print_solution(grid):
    for row in grid:
        print(','.join(row))

def is_valid_placement(grid, row, col, letter):
    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False
            
    return True

def solve(grid):
    # First, determine what letter should be on the minor diagonal
    diagonal_letters = set()
    for i in range(7):
        if grid[i][6-i] != '':
            diagonal_letters.add(grid[i][6-i])
    
    if len(diagonal_letters) > 1:
        return False  # Inconsistent diagonal
    
    diagonal_letter = list(diagonal_letters)[0] if diagonal_letters else 'c'  # Use 'c' if none found
    
    # Fill all empty diagonal positions with the diagonal letter
    for i in range(7):
        if grid[i][6-i] == '':
            if not is_valid_placement(grid, i, 6-i, diagonal_letter):
                return False
            grid[i][6-i] = diagonal_letter
    
    def solve_remaining(pos=0):
        if pos == 49:  # 7x7 = 49 cells
            return True
            
        row = pos // 7
        col = pos % 7
        
        # Skip filled cells and diagonal cells
        while pos < 49 and (grid[row][col] != '' or row + col == 6):
            pos += 1
            row = pos // 7
            col = pos % 7
        
        if pos == 49:
            return True
            
        for letter in 'abcdefg':
            if is_valid_placement(grid, row, col, letter):
                grid[row][col] = letter
                if solve_remaining(pos + 1):
                    return True
                grid[row][col] = ''
        return False
    
    return solve_remaining()

# Initial grid
grid = [
    ['c', 'g', '', 'd', 'e', 'b', ''],
    ['g', 'a', 'd', 'e', 'b', '', 'c'],
    ['', 'd', 'e', '', 'f', '', 'g'],
    ['d', 'e', '', 'f', '', '', ''],
    ['e', '', '', 'c', 'g', 'a', 'd'],
    ['b', 'f', 'c', 'g', 'a', '', ''],
    ['f', 'c', 'g', 'a', '', 'e', 'b']
]

if solve(grid):
    print_solution(grid)
else:
    print("No solution exists")