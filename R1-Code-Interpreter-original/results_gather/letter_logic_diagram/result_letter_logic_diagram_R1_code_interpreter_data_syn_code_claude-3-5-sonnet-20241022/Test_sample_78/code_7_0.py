def print_grid(grid):
    for row in grid:
        print(','.join(row))

def check_valid_placement(grid, row, col, letter):
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
    # Step 1: Verify and fill minor diagonal with 'd'
    for i in range(7):
        j = 6 - i
        if grid[i][j] != '' and grid[i][j] != 'd':
            return False
        grid[i][j] = 'd'
    
    # Step 2: Create list of empty positions (excluding minor diagonal)
    empty_positions = []
    for i in range(7):
        for j in range(7):
            if i + j != 6 and grid[i][j] == '':
                empty_positions.append((i, j))
    
    def backtrack(pos_index):
        if pos_index == len(empty_positions):
            return True
        
        row, col = empty_positions[pos_index]
        used_letters = set()
        
        # Get used letters in row and column
        for k in range(7):
            used_letters.add(grid[row][k])
            used_letters.add(grid[k][col])
        
        # Try each available letter
        available_letters = [l for l in 'abcdefg' if l not in used_letters]
        for letter in available_letters:
            if check_valid_placement(grid, row, col, letter):
                grid[row][col] = letter
                if backtrack(pos_index + 1):
                    return True
                grid[row][col] = ''
        
        return False
    
    return backtrack(0)

# Initialize grid
initial_grid = [
    ['', 'e', '', '', '', 'g', ''],
    ['e', 'c', 'a', 'b', '', 'd', 'f'],
    ['', '', 'b', 'g', '', 'f', 'e'],
    ['', 'b', 'g', '', '', 'e', 'c'],
    ['b', '', 'd', 'f', '', 'c', ''],
    ['g', 'd', '', '', 'c', 'a', 'b'],
    ['d', '', 'e', '', 'a', 'b', '']
]

# Create working copy
grid = [row[:] for row in initial_grid]

if solve(grid):
    print_grid(grid)
else:
    print("No solution found")