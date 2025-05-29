from itertools import permutations

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
    
    # Get remaining positions for each row
    for row in range(7):
        used = set(x for x in grid[row] if x != '')
        needed = set('abcdefg') - used
        empty_cols = [j for j in range(7) if grid[row][j] == '' and j != 6-row]  # exclude diagonal
        
        # Try all permutations of needed letters for empty positions
        for perm in permutations(needed, len(empty_cols)):
            # Check if this permutation would create valid columns
            valid = True
            for col, letter in zip(empty_cols, perm):
                # Check column validity
                col_letters = [grid[i][col] for i in range(7) if grid[i][col] != '']
                if letter in col_letters:
                    valid = False
                    break
            if valid:
                # Fill in the valid permutation
                for col, letter in zip(empty_cols, perm):
                    grid[row][j] = letter
                break
        else:
            return None
    
    return grid

result = solve()
if result:
    for row in result:
        print(','.join(row))
else:
    print("No solution exists")