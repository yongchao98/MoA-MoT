from itertools import permutations

def is_valid(grid, row, col, letter):
    # Check if the letter can be placed at grid[row][col]
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

def solve(grid, minor_diagonal_letter):
    # Try to fill the grid using permutations of 'abcdefg'
    for perm in permutations('abcdefg'):
        # Fill the grid row by row
        for i in range(7):
            for j in range(7):
                if grid[i][j] == '':
                    grid[i][j] = perm[j]
        
        # Check if the grid is valid
        valid = True
        for i in range(7):
            if not is_valid(grid, i, 0, grid[i][0]):
                valid = False
                break
        
        # Check the minor diagonal
        for i in range(7):
            if grid[i][6-i] != minor_diagonal_letter:
                valid = False
                break
        
        if valid:
            return True
        
        # Reset the grid
        for i in range(7):
            for j in range(7):
                if grid[i][j] not in 'abcdefg':
                    grid[i][j] = ''
    
    return False

def fill_grid():
    # Initial grid with given letters
    grid = [
        ['c', 'g', '', 'd', 'e', 'b', ''],
        ['g', 'a', 'd', 'e', 'b', '', 'c'],
        ['', 'd', 'e', '', 'f', '', 'g'],
        ['d', 'e', '', 'f', '', '', ''],
        ['e', '', '', 'c', 'g', 'a', 'd'],
        ['b', 'f', 'c', 'g', 'a', '', ''],
        ['f', 'c', 'g', 'a', '', 'e', 'b']
    ]

    # Determine the letter for the minor diagonal
    possible_letters = set('abcdefg')
    for i in range(7):
        for j in range(7):
            if i + j == 6 and grid[i][j] != '':
                possible_letters.intersection_update(grid[i][j])

    # Choose the letter for the minor diagonal
    minor_diagonal_letter = possible_letters.pop()

    # Solve the grid
    if solve(grid, minor_diagonal_letter):
        for row in grid:
            print(','.join(row))
    else:
        print("No solution found")

fill_grid()