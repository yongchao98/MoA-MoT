def print_grid(grid):
    for row in grid:
        print(row)

def get_prefilled_diagonal(grid):
    diagonal_cells = []
    for i in range(7):
        j = 6 - i
        if grid[i][j] != '':
            diagonal_cells.append(grid[i][j])
    return diagonal_cells

# Initial grid
grid = [
    ['', '', '', 'e', 'd', 'c', ''],
    ['f', '', 'e', '', '', '', ''],
    ['', '', 'd', '', 'g', '', ''],
    ['e', '', 'c', '', 'b', 'f', ''],
    ['d', '', '', 'b', '', 'a', 'e'],
    ['c', 'g', 'b', '', 'a', 'e', 'd'],
    ['', 'b', 'f', '', '', 'd', '']
]

# Let's analyze the pre-filled cells on the minor diagonal
diagonal_cells = get_prefilled_diagonal(grid)
print("Pre-filled cells on minor diagonal:", diagonal_cells)

# Looking at the output manually, I can see that 'g' is the only letter that can work
# for the minor diagonal given the constraints. Let me verify this:

def verify_g_diagonal():
    # Check if 'g' conflicts with any row or column where it would need to be placed
    # on the minor diagonal
    conflicts = []
    for i in range(7):
        j = 6 - i
        # Check if position already has a different letter
        if grid[i][j] != '' and grid[i][j] != 'g':
            conflicts.append(f"Position ({i},{j}) already has {grid[i][j]}")
        # Check row conflicts
        for k in range(7):
            if k != j and grid[i][k] == 'g':
                conflicts.append(f"Row {i} already has g at position {k}")
        # Check column conflicts
        for k in range(7):
            if k != i and grid[k][j] == 'g':
                conflicts.append(f"Column {j} already has g at position {k}")
    return conflicts

conflicts = verify_g_diagonal()
print("Conflicts for g on diagonal:", conflicts)