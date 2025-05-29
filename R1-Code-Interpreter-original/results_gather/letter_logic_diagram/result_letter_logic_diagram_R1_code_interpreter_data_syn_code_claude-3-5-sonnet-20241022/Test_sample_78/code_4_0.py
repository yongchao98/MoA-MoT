def print_grid(grid):
    for row in grid:
        print(','.join(row))

def get_missing_letters(used_letters):
    return [l for l in 'abcdefg' if l not in used_letters]

def is_valid_final(grid):
    # Check rows
    for row in grid:
        if sorted(row) != ['a','b','c','d','e','f','g']:
            return False
    
    # Check columns
    for j in range(7):
        col = [grid[i][j] for i in range(7)]
        if sorted(col) != ['a','b','c','d','e','f','g']:
            return False
    
    # Check minor diagonal
    for i in range(7):
        if grid[i][6-i] != 'd':
            return False
    
    return True

def solve(initial):
    # Step 1: Create a copy of initial grid
    grid = [row[:] for row in initial]
    
    # Step 2: Fill all minor diagonal cells with 'd'
    for i in range(7):
        if grid[i][6-i] == '':
            grid[i][6-i] = 'd'
        elif grid[i][6-i] != 'd':
            return None  # Invalid puzzle
    
    # Step 3: For each row, fill remaining cells
    for i in range(7):
        # Get used letters in this row
        used = [c for c in grid[i] if c != '']
        missing = get_missing_letters(used)
        
        # Fill empty cells with missing letters
        missing_idx = 0
        for j in range(7):
            if grid[i][j] == '':
                # Check if this letter would create a duplicate in column
                while missing_idx < len(missing):
                    letter = missing[missing_idx]
                    # Check if letter can be placed in this column
                    col_used = [grid[k][j] for k in range(7) if grid[k][j] != '']
                    if letter not in col_used:
                        grid[i][j] = letter
                        missing_idx += 1
                        break
                    missing_idx += 1
                if grid[i][j] == '':  # Couldn't place any letter
                    return None
    
    # Step 4: Verify solution
    if is_valid_final(grid):
        return grid
    return None

# Initialize grid with the given puzzle
initial_grid = [
    ['', 'e', '', '', '', 'g', 'd'],  # Force 'd' in top-right
    ['e', 'c', 'a', 'b', '', 'd', 'f'],
    ['', '', 'b', 'g', 'd', 'f', 'e'],
    ['', 'b', 'g', 'd', '', 'e', 'c'],
    ['b', '', 'd', 'f', '', 'c', ''],
    ['g', 'd', '', '', 'c', 'a', 'b'],
    ['d', '', 'e', '', 'a', 'b', '']
]

solution = solve(initial_grid)
if solution:
    print_grid(solution)
else:
    print("No solution found")