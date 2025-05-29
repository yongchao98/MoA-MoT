def is_valid_move(grid, row, col, letter, prefilled):
    # Check if violating prefilled position
    if prefilled[row][col] != '' and prefilled[row][col] != letter:
        return False

    # Check row
    if letter in grid[row]:
        return False

    # Check column
    if letter in [grid[i][col] for i in range(7)]:
        return False

    # Check minor diagonal
    if row + col == 6:  # if on minor diagonal
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '' and grid[i][j] != letter:
                return False

    return True

def find_next_empty(grid, prefilled):
    # First fill diagonal positions
    for i in range(7):
        j = 6 - i
        if grid[i][j] == '' and prefilled[i][j] == '':
            return (i, j)

    # Then fill other positions
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '' and prefilled[i][j] == '':
                return (i, j)
    return None

def solve_puzzle():
    prefilled = [
        ['', 'b', 'a', '', '', '', ''],
        ['b', '', '', 'd', '', '', ''],
        ['', '', 'd', 'f', '', '', 'b'],
        ['g', 'd', '', '', 'e', '', ''],
        ['d', 'f', 'c', 'e', 'b', 'a', ''],
        ['', '', '', '', '', '', 'd'],
        ['', 'e', '', '', '', '', 'f']
    ]

    # Initialize grid with prefilled positions
    grid = [[''] * 7 for _ in range(7)]
    for i in range(7):
        for j in range(7):
            if prefilled[i][j] != '':
                grid[i][j] = prefilled[i][j]

    def solve(grid):
        pos = find_next_empty(grid, prefilled)
        if not pos:
            return True

        row, col = pos

        # If on diagonal, try 'c' first
        if row + col == 6:
            # Check if any diagonal position is already filled
            diag_letter = None
            for i in range(7):
                j = 6 - i
                if grid[i][j] != '':
                    diag_letter = grid[i][j]
                    break
            
            if diag_letter:
                letters = [diag_letter]
            else:
                letters = ['c']  # Try 'c' first for diagonal
        else:
            letters = 'abcdefg'

        for letter in letters:
            if is_valid_move(grid, row, col, letter, prefilled):
                grid[row][col] = letter
                if solve(grid):
                    return True
                grid[row][col] = ''

        return False

    # Verify prefilled positions don't conflict
    for i in range(7):
        for j in range(7):
            if prefilled[i][j] != '':
                temp = grid[i][j]
                grid[i][j] = ''
                if not is_valid_move(grid, i, j, prefilled[i][j], prefilled):
                    print("Prefilled positions conflict!")
                    return None
                grid[i][j] = temp

    if solve(grid):
        return grid
    return None

# Solve and print
solution = solve_puzzle()
if solution:
    for row in solution:
        print(','.join(row))
else:
    print("No solution exists")