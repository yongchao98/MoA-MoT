def is_valid(grid, row, col, letter):
    # Check if the letter is already in the row
    if letter in grid[row]:
        return False
    # Check if the letter is already in the column
    for r in range(7):
        if grid[r][col] == letter:
            return False
    return True

def solve(grid, row=0, col=0):
    if row == 7:
        return True
    if col == 7:
        return solve(grid, row + 1, 0)
    if grid[row][col] != '':
        return solve(grid, row, col + 1)
    
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid, row, col + 1):
                return True
            grid[row][col] = ''
    return False

def fill_minor_diagonal(grid, letter):
    for i in range(7):
        grid[i][6-i] = letter

def main():
    grid = [
        ['','g','b','','','',''],
        ['','','','','c','',''],
        ['b','','','','','a',''],
        ['f','d','c','','','',''],
        ['','c','e','a','','b',''],
        ['','','','g','','f','d'],
        ['e','a','g','','f','','']
    ]
    
    for letter in 'abcdefg':
        # Check if the letter can be placed on the minor diagonal
        if all(is_valid(grid, r, 6-r, letter) for r in range(7)):
            fill_minor_diagonal(grid, letter)
            if solve(grid):
                for row in grid:
                    print(','.join(row))
                return
    print("No solution found")

main()