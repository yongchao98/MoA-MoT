def solve_puzzle(grid):
    # Helper function to find missing letters in a row or column
    def find_missing(letters):
        return set('abcdefg') - set(letters)

    # Check if placing a letter at grid[row][col] is valid
    def is_valid(row, col, letter):
        # Check row and column constraints
        if letter in grid[row] or letter in [grid[i][col] for i in range(7)]:
            return False
        return True

    # Backtracking function to fill the grid
    def backtrack(row, col):
        if row == 7:  # If we've filled all rows, the puzzle is solved
            return True
        if col == 7:  # Move to the next row
            return backtrack(row + 1, 0)
        if grid[row][col] != '':  # Skip pre-filled cells
            return backtrack(row, col + 1)

        # Try placing each letter in the current cell
        for letter in 'abcdefg':
            if is_valid(row, col, letter):
                grid[row][col] = letter
                if backtrack(row, col + 1):
                    return True
                grid[row][col] = ''  # Backtrack

        return False

    # Initialize the grid
    grid = [
        ['', 'a', 'g', '', '', '', ''],
        ['a', 'g', '', 'f', 'b', '', 'd'],
        ['', 'c', '', '', '', '', ''],
        ['', 'f', 'b', '', 'd', '', ''],
        ['', '', '', '', '', 'g', ''],
        ['e', '', '', '', 'g', 'c', ''],
        ['', '', 'a', 'g', 'c', '', 'b']
    ]

    # Choose a letter for the minor diagonal
    for letter in 'abcdefg':
        if all(is_valid(i, 6-i, letter) for i in range(7)):
            for i in range(7):
                grid[i][6-i] = letter
            break

    # Solve the puzzle using backtracking
    if backtrack(0, 0):
        # Format the output
        result = '\n'.join(','.join(row) for row in grid)
        return f"<<<\n{result}\n>>>"
    else:
        return "No solution found."

# Solve the puzzle and print the result
print(solve_puzzle([]))