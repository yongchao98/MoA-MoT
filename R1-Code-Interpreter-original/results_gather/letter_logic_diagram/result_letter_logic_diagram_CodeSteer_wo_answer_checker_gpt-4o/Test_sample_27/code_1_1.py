def solve_puzzle(grid):
    from itertools import permutations

    # Helper function to check if a letter can be placed at a given position
    def can_place(letter, row, col):
        return letter not in grid[row] and all(grid[r][col] != letter for r in range(7))

    # Try to fill the grid using backtracking
    def backtrack():
        for row in range(7):
            for col in range(7):
                if grid[row][col] == '':
                    for letter in 'abcdefg':
                        if can_place(letter, row, col):
                            grid[row][col] = letter
                            if backtrack():
                                return True
                            grid[row][col] = ''
                    return False
        return True

    # Fill the minor diagonal with the same letter
    for letter in 'abcdefg':
        if all(can_place(letter, r, 6 - r) for r in range(7)):
            for r in range(7):
                grid[r][6 - r] = letter
            if backtrack():
                break
            for r in range(7):
                grid[r][6 - r] = ''

    # Format the output
    return '\n'.join(','.join(row) for row in grid)

# Initial grid with given letters
initial_grid = [
    ['', 'd', 'f', 'e', 'b', 'a', 'g'],
    ['f', '', '', 'b', 'a', 'g', ''],
    ['e', 'b', '', '', '', 'c', 'd'],
    ['b', '', '', 'g', 'c', '', 'f'],
    ['b', 'a', '', 'c', 'd', '', ''],
    ['a', '', '', 'd', 'f', 'e', 'b'],
    ['', '', '', '', 'e', '', 'a']
]

# Solve the puzzle
solution = solve_puzzle(initial_grid)
print("<<<")
print(solution)
print(">>>")