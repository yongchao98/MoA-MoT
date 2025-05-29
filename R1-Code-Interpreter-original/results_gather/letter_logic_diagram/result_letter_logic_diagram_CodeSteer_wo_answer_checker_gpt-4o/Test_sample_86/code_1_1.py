def solve_puzzle(grid):
    from itertools import permutations

    def is_valid(grid):
        for i in range(7):
            if len(set(grid[i])) != 7 or len(set(row[i] for row in grid)) != 7:
                return False
        return True

    def fill_minor_diagonal(grid, letter):
        for i in range(7):
            grid[i][6-i] = letter

    def solve(grid, row=0):
        if row == 7:
            return is_valid(grid)
        
        missing_letters = set('abcdefg') - set(grid[row])
        for perm in permutations(missing_letters):
            new_row = list(grid[row])
            idx = 0
            for j in range(7):
                if new_row[j] == '':
                    new_row[j] = perm[idx]
                    idx += 1
            grid[row] = new_row
            if solve(grid, row + 1):
                return True
            grid[row] = list(grid[row])
        return False

    # Convert input grid to a list of lists
    grid = [list(row.split(',')) for row in grid.strip().split('\n')]

    # Try each letter for the minor diagonal
    for letter in 'abcdefg':
        temp_grid = [row[:] for row in grid]
        fill_minor_diagonal(temp_grid, letter)
        if solve(temp_grid):
            return '\n'.join(','.join(row) for row in temp_grid)

# Input grid
puzzle = """
,d,c,f,,b,a
,,f,e,b,a,
c,,e,b,a,,
,,b,a,g,d,
,,,g,d,c,
,,,d,,,e
,,,c,f,,
"""

# Solve the puzzle
solution = solve_puzzle(puzzle)
print(f"<<<\n{solution}\n>>>")