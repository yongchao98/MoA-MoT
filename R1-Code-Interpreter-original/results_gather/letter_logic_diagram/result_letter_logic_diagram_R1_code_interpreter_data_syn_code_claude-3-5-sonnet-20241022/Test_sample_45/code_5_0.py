def print_grid(grid):
    for row in grid:
        print(','.join(row))

def verify_solution(grid):
    # Verify rows
    for row in grid:
        if sorted(row) != ['a','b','c','d','e','f','g']:
            return False
    
    # Verify columns
    for j in range(7):
        col = [grid[i][j] for i in range(7)]
        if sorted(col) != ['a','b','c','d','e','f','g']:
            return False
    
    # Verify minor diagonal
    diag = [grid[i][6-i] for i in range(7)]
    if not all(x == 'b' for x in diag):
        return False
    
    return True

def matches_prefilled(grid, prefilled):
    for i in range(7):
        for j in range(7):
            if prefilled[i][j] != '' and grid[i][j] != prefilled[i][j]:
                return False
    return True

def solve():
    # Start with a valid Latin square and modify it to meet our constraints
    base = [
        ['c','d','a','e','f','b','g'],
        ['d','a','b','f','c','g','e'],
        ['a','e','f','b','g','c','d'],
        ['f','b','e','g','d','a','c'],
        ['e','f','g','c','b','d','a'],
        ['b','g','c','d','a','e','f'],
        ['g','c','d','a','e','f','b']
    ]
    
    # Force 'b' on minor diagonal
    for i in range(7):
        # Save current values
        j = 6-i
        temp = base[i][j]
        
        # Find where 'b' is in this row
        b_pos = base[i].index('b')
        
        # Swap 'b' to diagonal position
        base[i][j], base[i][b_pos] = base[i][b_pos], base[i][j]
        
        # Now need to fix the column that was affected
        # Find the row that has temp in the b_pos column
        for k in range(7):
            if k != i and base[k][b_pos] == temp:
                # Swap with the current value
                base[k][b_pos], base[i][b_pos] = base[i][b_pos], base[k][b_pos]
                break
    
    return base

# Initialize prefilled grid
prefilled = [
    ['c','','a','e','','b',''],
    ['','a','','f','','',''],
    ['','','f','b','g','',''],
    ['','','','','','',''],
    ['','','g','','','',''],
    ['b','','','d','a','e','f'],
    ['','','d','a','','f','b']
]

solution = solve()
print_grid(solution)