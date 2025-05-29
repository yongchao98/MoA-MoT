def print_grid(grid):
    for row in grid:
        print(','.join(row))

ORIGINAL = [
    ['a','e','','f','g','c',''],
    ['','b','','','','d','a'],
    ['b','f','g','c','d','a','e'],
    ['f','','c','','a','','b'],
    ['g','c','d','a','e','b','f'],
    ['c','','','e','b','','g'],
    ['d','','','','f','g','c']
]

def is_valid_row(row_num, row_values, grid):
    # Check if row matches original constraints
    for j in range(7):
        if ORIGINAL[row_num][j] != '' and row_values[j] != ORIGINAL[row_num][j]:
            return False
    
    # Check column uniqueness
    for j in range(7):
        letter = row_values[j]
        for i in range(row_num):
            if grid[i][j] == letter:
                return False
    
    # Check minor diagonal
    diag_letter = None
    for i in range(7):
        if i + (6-i) == 6 and grid[i][6-i] != '':
            diag_letter = grid[i][6-i]
            break
    
    if diag_letter:
        if row_num + (6-row_num) == 6 and row_values[6-row_num] != diag_letter:
            return False
    
    return True

def generate_valid_row(row_num, grid):
    from itertools import permutations
    
    # Get letters that must be in specific positions
    fixed_positions = [(j, ORIGINAL[row_num][j]) for j in range(7) if ORIGINAL[row_num][j] != '']
    
    # Get diagonal letter if it exists
    diag_letter = None
    for i in range(7):
        if i + (6-i) == 6 and grid[i][6-i] != '':
            diag_letter = grid[i][6-i]
            break
    
    if row_num + (6-row_num) == 6 and diag_letter:
        fixed_positions.append((6-row_num, diag_letter))
    
    # Try all permutations of remaining letters
    remaining_positions = [j for j in range(7) if not any(j == pos for pos, _ in fixed_positions)]
    remaining_letters = [l for l in 'abcdefg' if not any(l == letter for _, letter in fixed_positions)]
    
    for perm in permutations(remaining_letters):
        row = [''] * 7
        # Fill fixed positions
        for pos, letter in fixed_positions:
            row[pos] = letter
        # Fill remaining positions
        for pos, letter in zip(remaining_positions, perm):
            row[pos] = letter
        
        if is_valid_row(row_num, row, grid):
            return row
    
    return None

def solve():
    grid = []
    
    # Solve row by row
    for row_num in range(7):
        valid_row = generate_valid_row(row_num, grid)
        if valid_row is None:
            return False
        grid.append(valid_row)
    
    print_grid(grid)
    return True

solve()