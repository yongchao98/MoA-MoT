def verify_possible(grid):
    # Check if pre-filled positions allow for a consistent diagonal
    diagonal_letters = set()
    for i in range(7):
        if grid[i][6-i] != '':
            diagonal_letters.add(grid[i][6-i])
    
    # Print diagnostic information
    print("Pre-filled diagonal positions:")
    for i in range(7):
        if grid[i][6-i] != '':
            print(f"Position ({i},{6-i}): {grid[i][6-i]}")
    print(f"Unique letters on diagonal: {diagonal_letters}")
    
    # Check if any row/column already has multiple instances of same letter
    for i in range(7):
        row_letters = [x for x in grid[i] if x != '']
        if len(row_letters) != len(set(row_letters)):
            print(f"Row {i} has duplicate letters: {row_letters}")
            return False
    
    for j in range(7):
        col_letters = [grid[i][j] for i in range(7) if grid[i][j] != '']
        if len(col_letters) != len(set(col_letters)):
            print(f"Column {j} has duplicate letters: {col_letters}")
            return False
    
    return True

# Initialize grid
grid = [
    ['', 'a', '', 'g', 'b', '', ''],
    ['a', 'd', 'g', '', '', 'c', ''],
    ['d', '', '', 'e', 'c', '', 'a'],
    ['', 'b', 'e', '', 'f', 'a', 'd'],
    ['b', 'e', 'c', '', '', 'd', 'g'],
    ['', '', '', '', 'd', 'g', ''],
    ['', '', 'a', '', '', '', '']
]

print("Analyzing puzzle validity:")
verify_possible(grid)