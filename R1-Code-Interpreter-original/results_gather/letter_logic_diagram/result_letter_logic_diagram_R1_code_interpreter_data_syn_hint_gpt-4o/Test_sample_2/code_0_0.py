def solve_puzzle():
    # Initial grid setup
    grid = [
        ['','b','','f','g','','c'],
        ['b','e','','','a','c',''],
        ['','f','g','a','c','','b'],
        ['f','','a','','d','','e'],
        ['g','','','d','','e','f'],
        ['a','','','','','f','g'],
        ['c','d','','e','','g','a']
    ]
    
    # Determine the letter for the minor diagonal
    # Check which letter can fit in all diagonal positions
    possible_letters = set('abcdefg')
    for i in range(7):
        for j in range(7):
            if i + j == 6:  # Minor diagonal condition
                if grid[i][j] != '':
                    possible_letters.intersection_update(grid[i][j])
    
    # Choose the letter for the minor diagonal
    minor_diagonal_letter = possible_letters.pop()
    
    # Fill the minor diagonal
    for i in range(7):
        grid[i][6-i] = minor_diagonal_letter
    
    # Fill the rest of the grid
    for i in range(7):
        row_letters = set('abcdefg') - set(grid[i])
        for j in range(7):
            if grid[i][j] == '':
                col_letters = set(grid[k][j] for k in range(7))
                possible_letter = row_letters - col_letters
                grid[i][j] = possible_letter.pop()
    
    # Print the completed grid
    for row in grid:
        print(','.join(row))

solve_puzzle()