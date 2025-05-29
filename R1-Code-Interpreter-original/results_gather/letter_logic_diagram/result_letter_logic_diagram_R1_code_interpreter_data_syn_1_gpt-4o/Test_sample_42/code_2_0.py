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
    for letter in 'abcdefg':
        if all(can_place(letter, r, c) for r, c in minor_diagonal_positions):
            for r, c in minor_diagonal_positions:
                grid[r][c] = letter
            break

    # Backtracking function to fill the grid
    def backtrack(r, c):
        if r == 7:
            return True
        if c == 7:
            return backtrack(r + 1, 0)
        if grid[r][c] != '':
            return backtrack(r, c + 1)
        
        for letter in 'abcdefg':
            if can_place(letter, r, c):
                grid[r][c] = letter
                if backtrack(r, c + 1):
                    return True
                grid[r][c] = ''
        return False

    # Start backtracking from the first cell
    backtrack(0, 0)
    return grid

# Initial grid setup
grid = [
    ['', 'b', '', '', 'g', 'd', 'c'],
    ['', 'a', 'f', 'g', '', 'c', 'e'],
    ['a', 'f', '', '', 'c', 'e', 'b'],
    ['', 'g', 'd', '', 'e', '', ''],
    ['g', 'd', 'c', 'e', 'b', '', 'f'],
    ['', 'c', 'e', '', 'a', '', ''],
    ['', '', 'b', 'a', 'f', '', '']
]

# Solve the puzzle
solved_grid = solve_puzzle(grid)

# Format the output
for row in solved_grid:
    print(','.join(row))