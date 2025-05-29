def solve():
    # Initialize grid with givens
    grid = [
        ['', 'a', 'c', '', '', 'e', 'd'],
        ['a', '', 'f', 'g', '', 'd', 'b'],
        ['c', 'f', '', '', '', '', ''],
        ['', 'g', 'e', 'd', '', '', ''],
        ['', 'e', 'd', 'b', '', '', ''],
        ['', '', '', '', '', '', ''],
        ['d', '', 'a', '', '', '', '']
    ]
    
    # First, verify and fill minor diagonal with 'd'
    for i in range(7):
        if grid[i][6-i] != '' and grid[i][6-i] != 'd':
            return None
        grid[i][6-i] = 'd'
    
    def is_valid_column(col, letter, up_to_row):
        for i in range(up_to_row):
            if grid[i][col] == letter:
                return False
        return True
    
    def solve_row(row):
        if row == 7:
            return True
            
        # Get already used letters in this row
        used = set(x for x in grid[row] if x != '')
        needed = list(set('abcdefg') - used)
        empty_cols = [j for j in range(7) if grid[row][j] == '' and j != 6-row]  # exclude diagonal
        
        # Try all permutations of needed letters
        for perm in permutations(needed, len(empty_cols)):
            valid = True
            # Check if this permutation would create valid columns
            for col, letter in zip(empty_cols, perm):
                if not is_valid_column(col, letter, row):
                    valid = False
                    break
            
            if valid:
                # Fill in the valid permutation
                for col, letter in zip(empty_cols, perm):
                    grid[row][col] = letter
                if solve_row(row + 1):
                    return True
                # Backtrack
                for col in empty_cols:
                    grid[row][col] = ''
        return False
    
    if solve_row(0):
        return grid
    return None

result = solve()
if result:
    for row in result:
        print(','.join(row))
else:
    print("No solution exists")