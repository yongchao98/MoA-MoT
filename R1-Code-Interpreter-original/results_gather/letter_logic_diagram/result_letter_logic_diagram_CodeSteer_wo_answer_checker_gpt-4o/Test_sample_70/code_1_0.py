def construct_exact_cover_matrix(grid):
    size = 7
    letters = 'abcdefg'
    num_constraints = size * size * 3 + 7  # Rows, columns, and diagonal
    matrix = []

    # Map each cell to a unique constraint index
    def cell_constraint(row, col, letter):
        return row * size + col, size * size + row * size + (ord(letter) - ord('a')), 2 * size * size + col * size + (ord(letter) - ord('a'))

    # Add constraints for the minor diagonal
    def diagonal_constraint(row, col, letter):
        if row + col == size - 1:
            return 3 * size * size + (ord(letter) - ord('a'))
        return None

    # Create the exact cover matrix
    for row in range(size):
        for col in range(size):
            if grid[row][col] == '':
                for letter in letters:
                    row_constraints = [0] * num_constraints
                    cell_idx, row_idx, col_idx = cell_constraint(row, col, letter)
                    diag_idx = diagonal_constraint(row, col, letter)
                    row_constraints[cell_idx] = 1
                    row_constraints[row_idx] = 1
                    row_constraints[col_idx] = 1
                    if diag_idx is not None:
                        row_constraints[diag_idx] = 1
                    matrix.append(row_constraints)
                    # Debugging: Print the indices and row constraints
                    print(f"Row: {row}, Col: {col}, Letter: {letter}, Indices: {cell_idx}, {row_idx}, {col_idx}, {diag_idx}")
                    print(f"Row Constraints: {row_constraints}")

    return matrix

def solve_puzzle_with_dlx(grid):
    matrix = construct_exact_cover_matrix(grid)
    dlx = DLX(matrix)
    if dlx.search():
        # Convert the solution back to the grid format
        solution_grid = [['' for _ in range(7)] for _ in range(7)]
        for row in dlx.solution:
            r, c, l = row // 49, (row % 49) // 7, row % 7
            solution_grid[r][c] = 'abcdefg'[l]
        for row in solution_grid:
            print(','.join(row))

# Initial grid setup
grid = [
    ['', '', 'c', 'f', '', '', ''],
    ['e', 'c', 'f', '', '', '', ''],
    ['', 'f', 'g', 'b', 'a', '', 'e'],
    ['f', 'g', 'b', 'a', '', '', 'c'],
    ['g', 'b', 'a', 'd', 'e', '', 'f'],
    ['', '', '', '', 'c', '', ''],
    ['a', '', 'e', '', '', '', '']
]

solve_puzzle_with_dlx(grid)