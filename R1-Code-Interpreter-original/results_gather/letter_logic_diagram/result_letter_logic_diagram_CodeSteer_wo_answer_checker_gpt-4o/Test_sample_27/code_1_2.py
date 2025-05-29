def solve_puzzle(grid):
    # Initialize sets to track available letters for each row and column
    rows = [set('abcdefg') - set(row) for row in grid]
    cols = [set('abcdefg') - {grid[r][c] for r in range(7)} for c in range(7)]

    # Helper function to check if a letter can be placed at a given position
    def can_place(letter, row, col):
        return letter in rows[row] and letter in cols[col]

    # Try to fill the grid using backtracking
    def backtrack():
        for row in range(7):
            for col in range(7):
                if grid[row][col] == '':
                    for letter in rows[row] & cols[col]:
                        grid[row][col] = letter
                        rows[row].remove(letter)
                        cols[col].remove(letter)
                        if backtrack():
                            return True
                        grid[row][col] = ''
                        rows[row].add(letter)
                        cols[col].add(letter)
                    return False
        return True

    # Fill the minor diagonal with the same letter
    for letter in 'abcdefg':
        if all(can_place(letter, r, 6 - r) for r in range(7)):
            for r in range(7):
                grid[r][6 - r] = letter
                rows[r].remove(letter)
                cols[6 - r].remove(letter)
            if backtrack():
                break
            for r in range(7):
                grid[r][6 - r] = ''
                rows[r].add(letter)
                cols[6 - r].add(letter)

    # Format the output
    return '\n'.join(','.join(row) for row in grid)

# Initial grid with given letters
initial_grid = [
    ['', 'd', 'f', 'e', 'b', 'a', 'g'],
    ['f', '', '', 'b', 'a', 'g', ''],
    ['e', 'b', '', '', '', 'c', 'd'],
    ['b', '', '', 'g', 'c', '', 'f'],
    ['b', 'a', '', 'c', 'd', '', ''],
    ['a', '', '', 'd', 'f', 'e', 'b'],
    ['', '', '', '', 'e', '', 'a']
]

# Solve the puzzle
solution = solve_puzzle(initial_grid)
print("<<<")
print(solution)
print(">>>")