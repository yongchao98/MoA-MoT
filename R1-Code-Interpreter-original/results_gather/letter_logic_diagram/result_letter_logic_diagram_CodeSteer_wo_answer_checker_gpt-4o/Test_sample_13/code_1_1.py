def solve_puzzle(grid):
    from itertools import permutations

    # Helper function to check if a grid is valid
    def is_valid(grid):
        for i in range(7):
            if len(set(grid[i])) != 7:  # Check rows
                return False
            if len(set(grid[j][i] for j in range(7))) != 7:  # Check columns
                return False
        return True

    # Find the letter for the minor diagonal
    minor_diagonal = [grid[i][6-i] for i in range(7)]
    diagonal_letter = next((x for x in minor_diagonal if x), None)
    if not diagonal_letter:
        diagonal_letter = 'a'  # Default to 'a' if no pre-filled letter

    # Fill the minor diagonal with the determined letter
    for i in range(7):
        grid[i][6-i] = diagonal_letter

    # Backtracking function to fill the grid
    def backtrack(row):
        if row == 7:
            return is_valid(grid)
        if all(grid[row]):  # Skip filled rows
            return backtrack(row + 1)
        
        # Find missing letters in the current row
        missing = list(set('abcdefg') - set(grid[row]))
        for perm in permutations(missing):
            # Try filling the row with the permutation
            original_row = grid[row][:]
            perm_copy = list(perm)
            grid[row] = [grid[row][i] if grid[row][i] else perm_copy.pop(0) for i in range(7)]
            
            # Check column constraints before proceeding
            if all(len(set(grid[j][i] for j in range(row + 1))) == row + 1 for i in range(7)):
                if backtrack(row + 1):
                    return True
            
            grid[row] = original_row  # Backtrack
        return False

    # Start backtracking from the first row
    if backtrack(0):
        return '\n'.join(','.join(row) for row in grid)
    else:
        return "No solution found"

# Initial grid setup
grid = [
    ['a', '', 'd', '', 'e', 'c', ''],
    ['b', 'd', 'f', '', '', '', ''],
    ['d', '', '', 'c', '', '', 'a'],
    ['', '', 'c', 'g', '', '', ''],
    ['', '', 'g', '', 'b', '', ''],
    ['', 'g', '', 'b', '', '', 'f'],
    ['', '', '', 'd', '', 'e', '']
]

# Solve the puzzle
solution = solve_puzzle(grid)
print("<<<")
print(solution)
print(">>>")