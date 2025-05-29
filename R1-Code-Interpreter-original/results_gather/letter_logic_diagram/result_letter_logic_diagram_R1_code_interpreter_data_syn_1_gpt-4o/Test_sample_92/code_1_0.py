def solve_puzzle(grid):
    from collections import Counter
    
    # Function to check if a letter can be placed at a given position
    def can_place(letter, row, col):
        # Check row
        if letter in grid[row]:
            return False
        # Check column
        for r in range(7):
            if grid[r][col] == letter:
                return False
        return True

    # Fill the minor diagonal with a letter that can fit all positions
    minor_diagonal_positions = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    letters = set('abcdefg')
    
    # Try each letter for the minor diagonal
    for letter in letters:
        if all(can_place(letter, r, c) for r, c in minor_diagonal_positions):
            for r, c in minor_diagonal_positions:
                grid[r][c] = letter
            break

    # Fill the rest of the grid
    def backtrack():
        for r in range(7):
            for c in range(7):
                if grid[r][c] == '':
                    for letter in letters:
                        if can_place(letter, r, c):
                            grid[r][c] = letter
                            if backtrack():
                                return True
                            grid[r][c] = ''
                    return False
        return True

    backtrack()
    return grid

# Initial grid setup
grid = [
    ['a', 'f', 'b', 'd', '', 'c', ''],
    ['', '', 'd', '', 'c', '', 'a'],
    ['', '', 'e', '', '', '', 'f'],
    ['', '', 'c', '', '', 'f', ''],
    ['', '', '', 'a', '', 'b', ''],
    ['', '', '', 'f', 'b', '', ''],
    ['', 'a', 'f', '', '', 'e', '']
]

# Solve the puzzle
solved_grid = solve_puzzle(grid)

# Format the output
for row in solved_grid:
    print(','.join(row))