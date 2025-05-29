def solve_puzzle(grid):
    from collections import Counter
    
    # Function to find the missing letter in a row or column
    def find_missing(letters):
        all_letters = set('abcdefg')
        return list(all_letters - set(letters))
    
    # Backtracking function to fill the grid
    def backtrack(grid, row, col):
        if row == 7:
            return True  # All rows filled successfully
        
        if col == 7:
            return backtrack(grid, row + 1, 0)  # Move to the next row
        
        if grid[row][col] != '':
            return backtrack(grid, row, col + 1)  # Skip pre-filled cells
        
        # Find missing letters in the row and column
        row_missing = find_missing(grid[row])
        col_missing = find_missing([grid[x][col] for x in range(7)])
        
        # Find the intersection of missing letters in row and column
        possible_letters = set(row_missing) & set(col_missing)
        
        # If the cell is on the minor diagonal, it must match the diagonal letter
        if row + col == 6:
            possible_letters = {minor_diagonal_letter}
        
        for letter in possible_letters:
            grid[row][col] = letter
            if backtrack(grid, row, col + 1):
                return True
            grid[row][col] = ''  # Reset on backtrack
        
        return False
    
    # Determine the letter for the minor diagonal
    minor_diagonal_letter = 'd'  # Assume 'd' for now, adjust if needed
    
    # Start backtracking from the top-left corner
    if backtrack(grid, 0, 0):
        return grid
    else:
        return None

# Initial grid
grid = [
    ['', '', 'a', '', 'c', 'f', 'b'],
    ['g', 'a', '', '', '', 'b', 'e'],
    ['', '', 'c', 'f', 'b', 'e', 'g'],
    ['d', 'c', 'f', '', 'e', 'g', 'a'],
    ['', 'f', 'b', 'e', 'g', '', ''],
    ['f', 'b', 'e', 'g', 'a', 'd', ''],
    ['b', 'e', '', 'a', 'd', '', 'f']
]

# Solve the puzzle
solved_grid = solve_puzzle(grid)

# Print the solved grid
if solved_grid:
    for row in solved_grid:
        print(','.join(row))
else:
    print("No solution found")