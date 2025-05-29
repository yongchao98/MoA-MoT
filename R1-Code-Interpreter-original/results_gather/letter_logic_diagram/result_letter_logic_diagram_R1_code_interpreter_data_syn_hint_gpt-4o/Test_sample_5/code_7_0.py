def is_valid(grid, row, col, letter):
    # Check if the letter is already in the row
    if letter in grid[row]:
        return False
    # Check if the letter is already in the column
    for r in range(7):
        if grid[r][col] == letter:
            return False
    return True

def solve(grid, diagonal_letter):
    # Find the first empty cell
    for row in range(7):
        for col in range(7):
            if grid[row][col] == ' ':
                # Try each letter
                for letter in 'abcdefg':
                    if is_valid(grid, row, col, letter):
                        grid[row][col] = letter
                        if solve(grid, diagonal_letter):
                            return True
                        grid[row][col] = ' '
                return False
    # Check if the minor diagonal has the same letter
    for i in range(7):
        if grid[i][6-i] != diagonal_letter:
            return False
    return True

# Initial grid setup
grid = [
    [' ', 'f', 'b', 'e', 'g', 'd', 'a'],
    [' ', 'b', 'e', 'g', ' ', 'a', 'c'],
    [' ', ' ', ' ', 'd', 'a', ' ', ' '],
    ['e', ' ', ' ', 'a', 'c', 'f', ' '],
    ['g', ' ', 'a', 'c', 'f', 'b', ' '],
    ['d', 'a', 'c', 'f', 'b', 'e', 'g'],
    ['a', 'c', ' ', ' ', 'e', 'g', 'd']
]

# Try each letter for the minor diagonal
for diagonal_letter in 'abcdefg':
    # Place the letter on the minor diagonal
    for i in range(7):
        grid[i][6-i] = diagonal_letter
    
    # Solve the grid
    if solve(grid, diagonal_letter):
        for row in grid:
            print(','.join(row))
        break
else:
    print("No solution found")