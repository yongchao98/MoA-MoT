from copy import deepcopy

def solve_puzzle(grid):
    # Initialize available letters for each row and column
    row_available = [set('abcdefg') - set(row) for row in grid]
    col_available = [set('abcdefg') - set(grid[i][j] for i in range(7)) for j in range(7)]

    # Find the letter for the minor diagonal
    minor_diagonal_indices = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    used_letters = set(grid[i][j] for i, j in minor_diagonal_indices if grid[i][j] != '')
    all_letters = set('abcdefg')
    minor_diagonal_letter = (all_letters - used_letters).pop()

    # Fill the minor diagonal with the chosen letter
    for i, j in minor_diagonal_indices:
        grid[i][j] = minor_diagonal_letter
        row_available[i].discard(minor_diagonal_letter)
        col_available[j].discard(minor_diagonal_letter)

    # Try to fill the rest of the grid using backtracking with forward checking
    def backtrack():
        # Find the next empty cell
        for i in range(7):
            for j in range(7):
                if grid[i][j] == '':
                    # Try each available letter for this cell
                    for letter in row_available[i] & col_available[j]:
                        grid[i][j] = letter
                        row_available[i].remove(letter)
                        col_available[j].remove(letter)

                        if backtrack():
                            return True

                        # Backtrack
                        grid[i][j] = ''
                        row_available[i].add(letter)
                        col_available[j].add(letter)

                    return False
        return True

    # Start backtracking from the first cell
    backtrack()

    # Format the output
    return '\n'.join(','.join(row) for row in grid)

# Initial grid setup
grid = [
    ['', 'g', '', '', 'c', 'e', 'a'],
    ['g', '', 'b', '', '', '', ''],
    ['', '', 'c', 'e', '', '', ''],
    ['b', '', 'e', 'a', '', 'g', 'd'],
    ['', '', '', '', 'g', '', 'b'],
    ['', '', '', '', 'd', 'b', ''],
    ['', '', '', 'd', 'b', '', '']
]

# Solve the puzzle
solution = solve_puzzle(deepcopy(grid))
print("<<<")
print(solution)
print(">>>")