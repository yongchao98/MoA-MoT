def solve_puzzle(grid):
    from collections import defaultdict

    # Helper function to check if a letter can be placed at a given position
    def can_place(letter, row, col):
        return letter not in grid[row] and all(grid[r][col] != letter for r in range(7))

    # Find possible letters for the minor diagonal
    def find_possible_minor_diagonal_letters():
        possible_letters = set('abcdefg')
        for i in range(7):
            if grid[i][6 - i] != '':
                possible_letters.intersection_update(grid[i][6 - i])
        return possible_letters

    # Backtracking function with constraint propagation
    def backtrack(row=0, col=0):
        if row == 7:
            return True
        if col == 7:
            return backtrack(row + 1, 0)
        if grid[row][col] != '':
            return backtrack(row, col + 1)

        for letter in available_letters[row]:
            if can_place(letter, row, col):
                grid[row][col] = letter
                available_letters[row].remove(letter)
                if backtrack(row, col + 1):
                    return True
                grid[row][col] = ''
                available_letters[row].add(letter)
        return False

    # Initialize the grid
    grid = [list(row) for row in grid]

    # Determine possible minor diagonal letters
    possible_minor_diagonal_letters = find_possible_minor_diagonal_letters()

    # Try each possible minor diagonal letter
    for minor_diagonal_letter in possible_minor_diagonal_letters:
        # Set the minor diagonal
        for i in range(7):
            grid[i][6 - i] = minor_diagonal_letter

        # Track available letters for each row
        available_letters = [set('abcdefg') - set(row) for row in grid]

        # Solve the puzzle
        if backtrack():
            return '\n'.join(','.join(row) for row in grid)

    return "No solution found"

# Given puzzle grid
puzzle_grid = [
    ['b', 'c', '', '', '', '', 'f'],
    ['', '', '', '', 'a', '', 'b'],
    ['', '', 'd', '', 'f', 'b', 'c'],
    ['', 'd', '', '', 'b', 'c', ''],
    ['', 'a', 'f', '', '', '', 'g'],
    ['', 'f', '', 'c', 'e', '', ''],
    ['f', 'b', '', 'e', 'g', 'd', 'a']
]

# Solve the puzzle and print the solution
solution = solve_puzzle(puzzle_grid)
print("<<<")
print(solution)
print(">>>")