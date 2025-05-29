def solve_puzzle(grid):
    # Determine the letters present in each row and column
    letters = set('abcdefg')
    n = 7
    
    # Find the missing letters for each row and column
    row_missing = [letters - set(row) for row in grid]
    col_missing = [letters - set(grid[row][col] for row in range(n)) for col in range(n)]
    
    # Determine the letter for the minor diagonal
    # Check which letter can be placed in all diagonal positions
    for letter in letters:
        if all((grid[i][n-1-i] == '' or grid[i][n-1-i] == letter) for i in range(n)):
            minor_diagonal_letter = letter
            break
    
    # Fill the grid
    for i in range(n):
        for j in range(n):
            if grid[i][j] == '':
                if j == n-1-i:  # Minor diagonal
                    grid[i][j] = minor_diagonal_letter
                else:
                    # Find a letter that can be placed here
                    possible_letters = row_missing[i] & col_missing[j]
                    if len(possible_letters) == 1:
                        chosen_letter = possible_letters.pop()
                        grid[i][j] = chosen_letter
                        row_missing[i].remove(chosen_letter)
                        col_missing[j].remove(chosen_letter)
    
    return grid

# Initial grid with empty cells as ''
initial_grid = [
    ['', '', 'f', '', '', 'e', 'c'],
    ['d', 'f', 'b', 'g', 'e', 'c', ''],
    ['f', 'b', 'g', 'e', '', 'a', ''],
    ['', 'g', '', '', 'a', 'd', ''],
    ['', 'e', '', '', 'd', '', 'b'],
    ['e', '', 'a', 'd', 'f', 'b', 'g'],
    ['', '', 'd', 'f', '', 'g', '']
]

# Solve the puzzle
solved_grid = solve_puzzle(initial_grid)

# Print the solved grid
for row in solved_grid:
    print(','.join(row))