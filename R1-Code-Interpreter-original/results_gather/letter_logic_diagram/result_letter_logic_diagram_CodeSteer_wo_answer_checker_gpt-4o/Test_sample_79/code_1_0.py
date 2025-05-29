def is_valid(grid, row, col, letter):
    # Check if the letter can be placed at grid[row][col]
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

def find_minor_diagonal_letter(grid):
    # Count occurrences of each letter on the minor diagonal
    diagonal_counts = {}
    for i in range(7):
        letter = grid[i][6 - i]
        if letter:
            diagonal_counts[letter] = diagonal_counts.get(letter, 0) + 1
    # Choose the letter with the highest count
    return max(diagonal_counts, key=diagonal_counts.get)

def solve(grid, minor_diagonal_letter, row=0, col=0):
    if row == 7:
        return True  # Successfully filled the grid

    if col == 7:
        return solve(grid, minor_diagonal_letter, row + 1, 0)

    if grid[row][col] != '':
        return solve(grid, minor_diagonal_letter, row, col + 1)

    for letter in 'abcdefg':
        if (row + col == 6 and letter != minor_diagonal_letter) or not is_valid(grid, row, col, letter):
            continue

        grid[row][col] = letter
        if solve(grid, minor_diagonal_letter, row, col + 1):
            return True
        grid[row][col] = ''  # Backtrack

    return False

def solve_puzzle():
    grid = [
        ['d', 'e', 'c', 'a', 'f', '', 'g'],
        ['', 'c', '', '', 'b', 'g', ''],
        ['c', 'a', 'f', '', '', 'd', 'e'],
        ['a', '', '', 'g', 'd', 'e', 'c'],
        ['f', '', 'g', 'd', 'e', 'c', 'a'],
        ['b', 'g', '', '', 'c', '', 'f'],
        ['g', 'd', 'e', 'c', '', 'f', '']
    ]

    # Determine the letter for the minor diagonal
    minor_diagonal_letter = find_minor_diagonal_letter(grid)

    if solve(grid, minor_diagonal_letter):
        print("<<<")
        for row in grid:
            print(','.join(row))
        print(">>>")
    else:
        print("No solution found")

solve_puzzle()