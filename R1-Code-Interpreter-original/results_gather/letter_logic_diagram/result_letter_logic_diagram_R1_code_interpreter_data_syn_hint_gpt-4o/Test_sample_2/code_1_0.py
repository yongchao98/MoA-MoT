def is_valid(grid, row, col, letter):
    # Check if the letter can be placed at grid[row][col]
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

def solve(grid, minor_diagonal_letter, row=0, col=0):
    if row == 7:
        return True
    if col == 7:
        return solve(grid, minor_diagonal_letter, row + 1, 0)
    if grid[row][col] != '':
        return solve(grid, minor_diagonal_letter, row, col + 1)
    
    if row + col == 6:
        if is_valid(grid, row, col, minor_diagonal_letter):
            grid[row][col] = minor_diagonal_letter
            if solve(grid, minor_diagonal_letter, row, col + 1):
                return True
            grid[row][col] = ''
    else:
        for letter in 'abcdefg':
            if is_valid(grid, row, col, letter):
                grid[row][col] = letter
                if solve(grid, minor_diagonal_letter, row, col + 1):
                    return True
                grid[row][col] = ''
    return False

def solve_puzzle():
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
    possible_letters = set('abcdefg')
    for i in range(7):
        if grid[i][6-i] != '':
            possible_letters.intersection_update(grid[i][6-i])
    
    minor_diagonal_letter = possible_letters.pop()
    
    if solve(grid, minor_diagonal_letter):
        for row in grid:
            print(','.join(row))
    else:
        print("No solution found")

solve_puzzle()