def solve_puzzle(grid):
    # Determine the letter for the minor diagonal
    minor_diagonal_letter = None
    for i in range(7):
        if grid[i][6-i] != '':
            minor_diagonal_letter = grid[i][6-i]
            break
    
    # Fill the minor diagonal with the determined letter
    for i in range(7):
        grid[i][6-i] = minor_diagonal_letter
    
    # Function to find missing letters in a list
    def find_missing_letters(lst):
        return list(set('abcdefg') - set(lst))
    
    # Fill the grid
    for i in range(7):
        # Fill rows
        missing_in_row = find_missing_letters(grid[i])
        for j in range(7):
            if grid[i][j] == '':
                grid[i][j] = missing_in_row.pop()
    
    # Transpose the grid to work with columns
    transposed_grid = list(map(list, zip(*grid)))
    
    for j in range(7):
        # Fill columns
        missing_in_col = find_missing_letters(transposed_grid[j])
        for i in range(7):
            if transposed_grid[j][i] == '':
                transposed_grid[j][i] = missing_in_col.pop()
    
    # Transpose back to original grid
    solved_grid = list(map(list, zip(*transposed_grid)))
    
    # Format the output
    result = '\n'.join([','.join(row) for row in solved_grid])
    return result

# Initial grid
grid = [
    ['g', 'b', '', '', '', 'a', 'c'],
    ['b', 'e', 'f', 'd', 'a', 'c', 'g'],
    ['', '', '', '', 'c', 'g', ''],
    ['f', '', '', '', '', '', 'e'],
    ['', 'a', 'c', '', 'b', 'e', 'f'],
    ['a', '', 'g', 'b', 'e', 'f', ''],
    ['c', 'g', 'b', '', '', 'd', '']
]

# Solve the puzzle
print(solve_puzzle(grid))