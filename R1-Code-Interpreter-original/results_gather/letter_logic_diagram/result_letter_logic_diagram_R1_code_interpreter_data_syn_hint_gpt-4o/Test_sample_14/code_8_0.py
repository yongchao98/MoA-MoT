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
    
    # Function to check if a letter can be placed at grid[i][j]
    def can_place(i, j, letter):
        # Check row and column constraints
        for k in range(7):
            if grid[i][k] == letter or grid[k][j] == letter:
                return False
        return True
    
    # Backtracking function to fill the grid
    def backtrack(i, j):
        if i == 7:
            return True
        if j == 7:
            return backtrack(i + 1, 0)
        if grid[i][j] != '':
            return backtrack(i, j + 1)
        
        for letter in 'abcdefg':
            if can_place(i, j, letter):
                grid[i][j] = letter
                if backtrack(i, j + 1):
                    return True
                grid[i][j] = ''
        
        return False
    
    # Start backtracking from the first cell
    backtrack(0, 0)
    
    # Print the filled grid
    for row in grid:
        print(','.join(row))

fill_grid()