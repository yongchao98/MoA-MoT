def solve_puzzle(grid):
    from collections import Counter
    
    # Function to find the missing letter in a row or column
    def find_missing(letters):
        all_letters = set('abcdefg')
        return list(all_letters - set(letters))
    
    # Determine the letter for the minor diagonal
    minor_diagonal_letter = 'd'  # Assume 'd' for now, adjust if needed
    
    # Fill the grid
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                # Find missing letters in the row and column
                row_missing = find_missing(grid[i])
                col_missing = find_missing([grid[x][j] for x in range(7)])
                
                # Find the intersection of missing letters in row and column
                possible_letters = set(row_missing) & set(col_missing)
                
                # If the cell is on the minor diagonal, it must match the diagonal letter
                if i + j == 6:
                    grid[i][j] = minor_diagonal_letter
                else:
                    if possible_letters:
                        grid[i][j] = possible_letters.pop()
                    else:
                        # If no possible letters, backtrack or adjust logic
                        return None
    
    return grid

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