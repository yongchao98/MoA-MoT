def fill_grid():
    # Initial grid with given letters
    grid = [
        ['','f','','c','','',''],
        ['','d','','e','','',''],
        ['','','','','','b',''],
        ['c','e','g','','b','f',''],
        ['','','','','f','d',''],
        ['','a','b','','','c',''],
        ['','b','f','d','','','g']
    ]
    
    # Determine the letter for the minor diagonal
    # Check which letter is not present in any of the minor diagonal positions
    minor_diagonal_positions = [(0,6), (1,5), (2,4), (3,3), (4,2), (5,1), (6,0)]
    letters = set('abcdefg')
    
    # Find the letter for the minor diagonal
    for i, j in minor_diagonal_positions:
        if grid[i][j] in letters:
            letters.remove(grid[i][j])
    
    # The letter for the minor diagonal
    minor_diagonal_letter = letters.pop()
    
    # Fill the minor diagonal with the chosen letter
    for i, j in minor_diagonal_positions:
        grid[i][j] = minor_diagonal_letter
    
    # Function to find missing letters in a row
    def find_missing_letters(row):
        return list(set('abcdefg') - set(row))
    
    # Fill the grid
    for i in range(7):
        missing_letters = find_missing_letters(grid[i])
        for j in range(7):
            if grid[i][j] == '':
                # Find a letter that is not in the current column
                for letter in missing_letters:
                    if all(grid[k][j] != letter for k in range(7)):
                        grid[i][j] = letter
                        missing_letters.remove(letter)
                        break
    
    # Print the filled grid
    for row in grid:
        print(','.join(row))

fill_grid()