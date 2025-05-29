def is_valid(grid, row, col, letter):
    # Check if the letter can be placed at grid[row][col]
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

def solve(grid, row=0, col=0, minor_diagonal_letter=None):
    if row == 7:
        # If we've filled all rows, return the solution
        return grid
    
    if col == 7:
        # Move to the next row
        return solve(grid, row + 1, 0, minor_diagonal_letter)
    
    if grid[row][col] != '':
        # Skip pre-filled cells
        return solve(grid, row, col + 1, minor_diagonal_letter)
    
    if row + col == 6:
        # Fill the minor diagonal with the chosen letter
        if minor_diagonal_letter is None:
            # Try each letter for the minor diagonal
            for letter in 'abcdefg':
                if is_valid(grid, row, col, letter):
                    grid[row][col] = letter
                    result = solve(grid, row, col + 1, letter)
                    if result:
                        return result
                    grid[row][col] = ''
            return None
        else:
            if is_valid(grid, row, col, minor_diagonal_letter):
                grid[row][col] = minor_diagonal_letter
                result = solve(grid, row, col + 1, minor_diagonal_letter)
                if result:
                    return result
                grid[row][col] = ''
            return None
    
    # Try each letter for the current cell
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            result = solve(grid, row, col + 1, minor_diagonal_letter)
            if result:
                return result
            grid[row][col] = ''
    
    return None

def solve_puzzle(puzzle):
    # Convert the grid into a list of lists for easier manipulation
    grid = [row.split(',') for row in puzzle.strip().split('\n')]
    
    # Solve the puzzle
    solution = solve(grid)
    
    # Format the output
    if solution:
        result = '\n'.join(','.join(row) for row in solution)
        return f"<<<\n{result}\n>>>"
    else:
        return "No solution found."

# Initial grid
puzzle = """
b,,,c,g,f,a
,e,c,g,f,a,b
e,c,g,f,a,b,d
c,,f,,b,,e
g,f,,b,d,,
,,b,,e,c,g
,,d,e,,g,f
"""

# Solve the puzzle
print(solve_puzzle(puzzle))