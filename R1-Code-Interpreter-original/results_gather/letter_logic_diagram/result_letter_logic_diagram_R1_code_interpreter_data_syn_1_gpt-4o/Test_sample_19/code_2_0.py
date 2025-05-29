def solve_puzzle():
    from itertools import permutations

    # Initialize the grid
    grid = [
        ['g', 'b', 'f', '', 'c', 'a', ''],
        ['', 'f', 'd', '', 'a', 'e', ''],
        ['', 'd', 'c', '', '', '', ''],
        ['', '', '', '', '', '', ''],
        ['', '', 'e', '', 'b', '', ''],
        ['a', 'e', '', '', '', 'd', 'c'],
        ['e', '', '', '', '', 'c', '']
    ]

    # Function to check if the grid is valid
    def is_valid(grid):
        for i in range(7):
            if len(set(grid[i])) != 7 or len(set(row[i] for row in grid)) != 7:
                return False
        return True

    # Try all permutations for the minor diagonal
    for letter in 'abcdefg':
        # Fill the minor diagonal
        for i in range(7):
            grid[i][6-i] = letter

        # Try to fill the rest of the grid
        for perm in permutations('abcdefg'):
            for i in range(7):
                if grid[i].count('') > 0:
                    grid[i] = [perm[j] if cell == '' else cell for j, cell in enumerate(grid[i])]
            if is_valid(grid):
                break

    # Format the output
    result = '\n'.join([','.join(row) for row in grid])
    return result

# Solve the puzzle
solution = solve_puzzle()
print(solution)