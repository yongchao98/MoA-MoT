def solve_puzzle(grid):
    def is_valid(grid, row, col, letter):
        # Check if the letter is not in the current row and column
        for i in range(7):
            if grid[row][i] == letter or grid[i][col] == letter:
                return False
        # Check minor diagonal constraint
        if row + col == 6:
            for i in range(7):
                if grid[i][6-i] != '' and grid[i][6-i] != letter:
                    return False
        return True

    def solve(grid):
        for row in range(7):
            for col in range(7):
                if grid[row][col] == '':
                    for letter in 'abcdefg':
                        if is_valid(grid, row, col, letter):
                            grid[row][col] = letter
                            if solve(grid):
                                return True
                            grid[row][col] = ''  # Backtrack
                    return False
        return True

    # Parse the input grid
    grid = [row.split(',') for row in grid.strip().split('\n')]
    
    # Solve the puzzle
    solve(grid)
    
    # Format the output
    result = '\n'.join([','.join(row) for row in grid])
    return f"<<<\n{result}\n>>>"

# Input grid
input_grid = """
c,,g,,f,e,
a,g,,,,d,
,,,,d,,
,,,d,,a,g
,e,,c,a,g,b
e,,,a,g,b,f
d,,a,,b,,
"""

# Solve and print the result
print(solve_puzzle(input_grid))