def solve_puzzle(grid):
    from collections import defaultdict

    def get_possibilities(grid):
        possibilities = defaultdict(set)
        for i in range(7):
            for j in range(7):
                if grid[i][j] == '':
                    row_letters = set(grid[i])
                    col_letters = set(row[j] for row in grid)
                    possibilities[(i, j)] = set('abcdefg') - row_letters - col_letters
        return possibilities

    def select_minor_diagonal_letter(grid):
        diagonal_letters = [grid[i][6-i] for i in range(7) if grid[i][6-i] != '']
        if diagonal_letters:
            return diagonal_letters[0]
        return 'a'  # Default if no letters are pre-filled

    def is_valid(grid):
        for i in range(7):
            if len(set(grid[i])) != 7 or len(set(row[i] for row in grid)) != 7:
                return False
        minor_diagonal_letter = grid[0][6]
        for i in range(7):
            if grid[i][6-i] != minor_diagonal_letter:
                return False
        return True

    def fill_grid(grid, possibilities):
        if all(grid[i][j] != '' for i in range(7) for j in range(7)):
            return is_valid(grid)

        # Find the cell with the fewest possibilities
        min_pos = None
        min_count = float('inf')
        for pos, opts in possibilities.items():
            if len(opts) < min_count:
                min_count = len(opts)
                min_pos = pos

        if min_pos is None:
            return False

        i, j = min_pos
        for letter in possibilities[min_pos]:
            new_grid = [list(row) for row in grid]
            new_grid[i][j] = letter
            new_possibilities = get_possibilities(new_grid)
            if fill_grid(new_grid, new_possibilities):
                for r in range(7):
                    grid[r] = new_grid[r]
                return True
        return False

    # Select the minor diagonal letter
    minor_diagonal_letter = select_minor_diagonal_letter(grid)
    for i in range(7):
        if grid[i][6-i] == '':
            grid[i][6-i] = minor_diagonal_letter

    possibilities = get_possibilities(grid)
    if not fill_grid(grid, possibilities):
        raise ValueError("No solution found")

    return grid

# Initial grid with empty cells as ''
grid = [
    ['c', 'b', '', '', 'd', '', ''],
    ['b', 'g', 'f', '', 'e', '', 'c'],
    ['g', '', 'd', 'e', 'a', 'c', 'b'],
    ['f', 'd', 'e', '', 'c', 'b', 'g'],
    ['d', 'e', '', '', 'b', 'g', ''],
    ['e', 'a', 'c', 'b', 'g', 'f', 'd'],
    ['', '', 'b', 'g', '', '', 'e']
]

solution = solve_puzzle(grid)
for row in solution:
    print(','.join(row))