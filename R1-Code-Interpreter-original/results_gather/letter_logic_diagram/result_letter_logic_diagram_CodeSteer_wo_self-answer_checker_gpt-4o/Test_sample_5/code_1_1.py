def solve_puzzle(grid):
    from copy import deepcopy

    # Helper function to check if the grid is valid
    def is_valid(grid):
        for i in range(7):
            if len(set(grid[i])) != 7:  # Check rows
                return False
            if len(set(grid[j][i] for j in range(7))) != 7:  # Check columns
                return False
        return True

    # Find the letter for the minor diagonal
    diagonal_indices = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    diagonal_letters = set(grid[i][j] for i, j in diagonal_indices if grid[i][j] != '')
    all_letters = set('abcdefg')
    diagonal_letter = (all_letters - diagonal_letters).pop()

    # Fill the diagonal with the chosen letter
    for i, j in diagonal_indices:
        grid[i][j] = diagonal_letter

    # Function to get possible letters for a cell
    def get_possible_letters(grid, row, col):
        if grid[row][col] != '':
            return set()
        row_letters = set(grid[row])
        col_letters = set(grid[i][col] for i in range(7))
        return all_letters - row_letters - col_letters

    # Try to fill the grid using backtracking with forward checking
    def backtrack(grid, row=0):
        if row == 7:
            return is_valid(grid)
        for col in range(7):
            if grid[row][col] == '':
                possible_letters = get_possible_letters(grid, row, col)
                for letter in possible_letters:
                    grid[row][col] = letter
                    if backtrack(grid, row):
                        return True
                    grid[row][col] = ''  # Reset cell
                return False
        return backtrack(grid, row + 1)

    # Start backtracking from the first row
    if backtrack(grid):
        # Print the solution
        for row in grid:
            print(','.join(row))
    else:
        print("No solution found")

# Initial grid with empty cells as ''
grid = [
    ['', 'f', 'b', 'e', 'g', 'd', 'a'],
    ['b', 'e', 'g', '', 'a', 'c', ''],
    ['', '', '', 'd', 'a', '', ''],
    ['e', '', '', 'a', 'c', 'f', ''],
    ['g', '', 'a', 'c', 'f', 'b', ''],
    ['d', 'a', 'c', 'f', 'b', 'e', 'g'],
    ['a', 'c', '', '', 'e', 'g', 'd']
]

solve_puzzle(grid)