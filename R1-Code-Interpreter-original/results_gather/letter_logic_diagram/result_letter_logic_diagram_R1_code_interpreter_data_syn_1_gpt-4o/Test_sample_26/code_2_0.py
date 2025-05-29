def solve_puzzle(grid):
    from collections import Counter

    # Initialize the grid
    n = 7
    letters = set('abcdefg')
    minor_diagonal_indices = [(i, n - 1 - i) for i in range(n)]

    # Determine the missing letters for each row and column
    row_missing = [letters - set(row) for row in grid]
    col_missing = [letters - set(grid[i][j] for i in range(n)) for j in range(n)]

    # Determine the letter for the minor diagonal
    diagonal_letters = [grid[i][n - 1 - i] for i in range(n) if grid[i][n - 1 - i] != '']
    diagonal_letter = Counter(diagonal_letters).most_common(1)[0][0] if diagonal_letters else 'a'

    # Fill the grid
    for i in range(n):
        for j in range(n):
            if grid[i][j] == '':
                if (i, j) in minor_diagonal_indices:
                    grid[i][j] = diagonal_letter
                else:
                    # Find a letter that can be placed here
                    possible_letters = row_missing[i] & col_missing[j]
                    if possible_letters:
                        chosen_letter = possible_letters.pop()
                        grid[i][j] = chosen_letter
                        row_missing[i].remove(chosen_letter)
                        col_missing[j].remove(chosen_letter)

    # Format the output
    return '\n'.join(','.join(row) for row in grid)

# Input grid
grid = [
    ['', 'g', 'a', '', 'd', '', ''],
    ['a', '', '', '', 'c', 'f', ''],
    ['a', '', 'd', 'b', 'c', '', 'g'],
    ['e', 'd', '', 'c', 'f', 'g', 'a'],
    ['', 'b', 'c', '', '', 'a', ''],
    ['', 'c', '', 'g', 'a', 'e', 'd'],
    ['', 'f', '', 'a', '', 'd', 'b']
]

# Solve the puzzle
solution = solve_puzzle(grid)
print(solution)