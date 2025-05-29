def solve_puzzle(grid):
    n = 7
    letters = set('abcdefg')
    
    # Determine the missing letters for each row and column
    row_missing = [letters - set(row) for row in grid]
    col_missing = [letters - set(grid[row][col] for row in range(n)) for col in range(n)]
    
    # Find the letter for the minor diagonal
    # Check which letter can be placed on the diagonal
    for letter in letters:
        if all(letter in row_missing[i] or grid[i][n-1-i] == letter for i in range(n)):
            diagonal_letter = letter
            break
    
    # Fill the grid
    for i in range(n):
        for j in range(n):
            if grid[i][j] == '':
                if j == n-1-i:  # If it's on the minor diagonal
                    grid[i][j] = diagonal_letter
                else:
                    # Find a letter that can be placed here
                    possible_letters = row_missing[i] & col_missing[j]
                    grid[i][j] = possible_letters.pop()
                    row_missing[i].remove(grid[i][j])
                    col_missing[j].remove(grid[i][j])
    
    return grid

# Initial grid
grid = [
    ['e', 'a', '', '', '', 'c', 'f'],
    ['', 'b', '', '', 'c', 'f', ''],
    ['b', 'g', 'd', 'c', 'f', '', 'a'],
    ['g', 'd', '', 'f', 'e', '', 'b'],
    ['d', '', '', 'e', '', '', 'g'],
    ['', 'f', 'e', '', 'b', '', 'd'],
    ['', '', 'a', '', 'g', 'd', 'c']
]

# Solve the puzzle
solved_grid = solve_puzzle(grid)

# Print the solved grid
for row in solved_grid:
    print(','.join(row))