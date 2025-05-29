def solve_puzzle(grid):
    # Helper function to find missing letters in a list
    def find_missing(letters):
        return [chr(i) for i in range(ord('a'), ord('g') + 1) if chr(i) not in letters]

    # Initialize sets to track available letters for rows, columns, and diagonal
    available_rows = [set(find_missing(row)) for row in grid]
    available_cols = [set(find_missing([grid[i][j] for i in range(7)])) for j in range(7)]
    diagonal_positions = [(i, 6 - i) for i in range(7)]
    
    # Determine a valid diagonal letter
    diagonal_letter = None
    for letter in 'abcdefg':
        if all(grid[i][6 - i] in ('', letter) for i in range(7)):
            diagonal_letter = letter
            break

    # Backtracking function
    def backtrack(row, col):
        if row == 7:
            return True
        if col == 7:
            return backtrack(row + 1, 0)
        if grid[row][col] != '':
            return backtrack(row, col + 1)

        for letter in 'abcdefg':
            if (letter in available_rows[row] and
                letter in available_cols[col] and
                (row, col) not in diagonal_positions or letter == diagonal_letter):
                
                # Place the letter
                grid[row][col] = letter
                available_rows[row].remove(letter)
                available_cols[col].remove(letter)

                # Move to the next cell
                if backtrack(row, col + 1):
                    return True

                # Backtrack
                grid[row][col] = ''
                available_rows[row].add(letter)
                available_cols[col].add(letter)

        return False

    # Start backtracking from the first cell
    backtrack(0, 0)

    # Format the output
    result = '\n'.join([','.join(row) for row in grid])
    return f"<<<\n{result}\n>>>"

# Initial grid
grid = [
    ['e', 'f', '', 'd', '', 'a', 'c'],
    ['', 'b', '', '', 'a', '', 'e'],
    ['', 'd', 'g', 'a', '', 'e', ''],
    ['d', 'g', 'a', 'c', 'e', 'f', 'b'],
    ['g', 'a', '', 'e', 'f', 'b', 'd'],
    ['', 'c', '', '', '', 'd', 'g'],
    ['', 'e', 'f', 'b', 'd', 'g', 'a']
]

# Solve the puzzle
print(solve_puzzle(grid))