def solve():
    # Initialize grid
    grid = [['' for _ in range(7)] for _ in range(7)]
    
    # Fill in givens
    givens = [
        (0,1,'a'), (0,2,'c'), (0,5,'e'), (0,6,'d'),
        (1,0,'a'), (1,2,'f'), (1,3,'g'), (1,5,'d'), (1,6,'b'),
        (2,0,'c'), (2,1,'f'),
        (3,1,'g'), (3,2,'e'), (3,3,'d'),
        (4,1,'e'), (4,2,'d'), (4,3,'b'),
        (6,0,'d'), (6,2,'a')
    ]
    
    # Fill givens and verify they don't conflict
    for r,c,v in givens:
        if grid[r][c] != '':
            return None
        grid[r][c] = v
    
    # FORCE all minor diagonal positions to be 'd'
    for i in range(7):
        if grid[i][6-i] != '' and grid[i][6-i] != 'd':
            return None
        grid[i][6-i] = 'd'
    
    def is_safe(r, c, v):
        if r + c == 6:  # Minor diagonal
            return v == 'd'
        
        # Check row
        if v in grid[r]:
            return False
            
        # Check column
        if v in [grid[i][c] for i in range(7)]:
            return False
            
        return True
    
    def find_empty():
        # First ensure minor diagonal is filled
        for i in range(7):
            if grid[i][6-i] == '':
                return (i, 6-i)
        # Then fill other cells
        for r in range(7):
            for c in range(7):
                if grid[r][c] == '':
                    return (r, c)
        return None
    
    def backtrack():
        pos = find_empty()
        if not pos:
            return True
            
        r, c = pos
        candidates = ['d'] if r + c == 6 else 'abcdefg'
        
        for v in candidates:
            if is_safe(r, c, v):
                grid[r][c] = v
                if backtrack():
                    return True
                grid[r][c] = ''
        return False
    
    if backtrack():
        return grid
    return None

result = solve()
if result:
    for row in result:
        print(','.join(row))
else:
    print("No solution exists")