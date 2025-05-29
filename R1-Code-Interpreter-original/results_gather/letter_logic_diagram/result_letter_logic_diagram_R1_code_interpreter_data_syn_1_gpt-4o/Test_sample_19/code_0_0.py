def solve_puzzle(grid):
    from collections import Counter

    # Function to find missing letters in a row or column
    def find_missing(letters):
        all_letters = set('abcdefg')
        return list(all_letters - set(letters))

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

    # Determine the missing letters for each row
    row_missing = [find_missing(row) for row in grid]

    # Determine the missing letters for each column
    columns = [[grid[row][col] for row in range(7)] for col in range(7)]
    col_missing = [find_missing(col) for col in columns]

    # Try to fill the minor diagonal with a consistent letter
    for letter in 'abcdefg':
        if all(letter in row_missing[i] or grid[i][6-i] == letter for i in range(7)):
            # Fill the minor diagonal
            for i in range(7):
                if grid[i][6-i] == '':
                    grid[i][6-i] = letter
                    row_missing[i].remove(letter)
                    col_missing[6-i].remove(letter)
            break

    # Fill the remaining cells
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                # Find a letter that can be placed here
                for letter in row_missing[i]:
                    if letter in col_missing[j]:
                        grid[i][j] = letter
                        row_missing[i].remove(letter)
                        col_missing[j].remove(letter)
                        break

    # Format the output
    result = '\n'.join([','.join(row) for row in grid])
    return result

# Solve the puzzle
solution = solve_puzzle([])
print(solution)