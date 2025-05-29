def verify_constraints(grid, row, col, val, diag_letter):
    # Check row
    if val in grid[row]:
        return False
    
    # Check column
    for i in range(7):
        if grid[i][col] == val:
            return False
    
    # Check minor diagonal
    if row + col == 6 and val != diag_letter:
        return False
        
    return True

def solve_with_diagonal(grid, diag_letter, filled_positions):
    if len(filled_positions) == 49:  # All cells filled
        return True
    
    # Get next unfilled position
    pos = len(filled_positions)
    row, col = pos // 7, pos % 7
    
    # If this is a diagonal position, must be diag_letter
    if row + col == 6:
        if verify_constraints(grid, row, col, diag_letter, diag_letter):
            grid[row][col] = diag_letter
            filled_positions.append((row, col))
            if solve_with_diagonal(grid, diag_letter, filled_positions):
                return True
            grid[row][col] = ''
            filled_positions.pop()
    else:
        for val in 'abcdefg':
            if verify_constraints(grid, row, col, val, diag_letter):
                grid[row][col] = val
                filled_positions.append((row, col))
                if solve_with_diagonal(grid, diag_letter, filled_positions):
                    return True
                grid[row][col] = ''
                filled_positions.pop()
    
    return False

# Create empty grid
grid = [['' for _ in range(7)] for _ in range(7)]

# Pre-fill required positions
prefilled = [
    (0,0,'a'), (0,1,'g'), (0,4,'f'), (0,6,'b'),
    (1,1,'d'), (1,4,'e'), (1,6,'a'),
    (2,0,'d'), (2,2,'f'), (2,3,'e'), (2,5,'a'), (2,6,'g'),
    (3,0,'c'), (3,1,'f'), (3,2,'e'), (3,3,'b'), (3,4,'a'), (3,5,'g'),
    (4,2,'b'), (4,4,'g'), (4,5,'d'), (4,6,'c'),
    (5,0,'e'), (5,1,'b'), (5,2,'a'), (5,3,'g'), (5,4,'d'), (5,5,'c'),
    (6,0,'b'), (6,2,'g'), (6,3,'d'), (6,4,'c'), (6,5,'f'), (6,6,'e')
]

# Try each possible diagonal letter
for diag_letter in 'abcdefg':
    # Reset grid
    grid = [['' for _ in range(7)] for _ in range(7)]
    
    # Verify pre-filled positions with this diagonal letter
    valid = True
    for row, col, val in prefilled:
        if row + col == 6 and val != diag_letter:
            valid = False
            break
        if not verify_constraints(grid, row, col, val, diag_letter):
            valid = False
            break
        grid[row][col] = val

    if valid:
        filled_positions = [(r,c) for r,c,_ in prefilled]
        filled_positions.sort(key=lambda x: x[0]*7 + x[1])
        
        if solve_with_diagonal(grid, diag_letter, filled_positions):
            print('<<<')
            for row in grid:
                print(','.join(row))
            print('>>>')
            break