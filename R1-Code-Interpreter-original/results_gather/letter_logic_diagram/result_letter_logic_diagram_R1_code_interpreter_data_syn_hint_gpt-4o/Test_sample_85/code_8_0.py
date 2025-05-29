from collections import Counter

def is_valid(grid, row, col, letter):
    # Check if the letter can be placed at grid[row][col]
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

def solve(grid, empty_cells, index, diagonal_letter):
    if index == len(empty_cells):
        return True  # All cells are filled

    row, col = empty_cells[index]

    # If the cell is on the minor diagonal, it must be the diagonal letter
    if (row, col) in [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]:
        if is_valid(grid, row, col, diagonal_letter):
            grid[row][col] = diagonal_letter
            if solve(grid, empty_cells, index + 1, diagonal_letter):
                return True
            grid[row][col] = ''  # Backtrack
    else:
        for letter in 'abcdefg':
            if is_valid(grid, row, col, letter):
                grid[row][col] = letter
                if solve(grid, empty_cells, index + 1, diagonal_letter):
                    return True
                grid[row][col] = ''  # Backtrack

    return False

def solve_puzzle(grid):
    # Identify empty cells
    empty_cells = [(i, j) for i in range(7) for j in range(7) if grid[i][j] == '']

    # Determine the letter for the minor diagonal
    # Count occurrences of each letter in the minor diagonal
    diagonal_letters = [grid[i][j] for i, j in [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)] if grid[i][j] != '']
    letter_counts = Counter(diagonal_letters)
    
    # Choose a letter for the diagonal that is not over-represented
    for letter in 'abcdefg':
        if letter_counts[letter] == 0:
            diagonal_letter = letter
            break

    # Solve the puzzle using backtracking
    if solve(grid, empty_cells, 0, diagonal_letter):
        return '\n'.join(','.join(row) for row in grid)
    else:
        return "No solution found"

# Initial grid
grid = [
    ['b', 'c', '', '', '', '', 'f'],
    ['', '', '', '', 'a', '', 'b'],
    ['', '', 'd', '', 'f', 'b', 'c'],
    ['', 'd', '', '', 'b', 'c', ''],
    ['', 'a', 'f', '', '', '', 'g'],
    ['', 'f', '', 'c', 'e', '', ''],
    ['f', 'b', '', 'e', 'g', 'd', 'a']
]

# Solve the puzzle
solution = solve_puzzle(grid)
print(solution)