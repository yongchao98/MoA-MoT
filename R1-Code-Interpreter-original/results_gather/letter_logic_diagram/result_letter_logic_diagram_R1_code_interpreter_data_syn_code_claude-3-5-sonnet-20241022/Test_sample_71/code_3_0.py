def print_solution(grid):
    for row in grid:
        print(','.join(row))

def is_valid_complete_solution(grid):
    # Check rows
    for row in grid:
        if sorted(row) != ['a','b','c','d','e','f','g']:
            return False
    
    # Check columns
    for j in range(7):
        col = [grid[i][j] for i in range(7)]
        if sorted(col) != ['a','b','c','d','e','f','g']:
            return False
    
    # Check minor diagonal
    minor_diag = [grid[i][6-i] for i in range(7)]
    if len(set(minor_diag)) != 1:  # All elements should be same
        return False
    
    # Check if solution matches initial constraints
    initial = [
        ['', 'd', '', '', 'b', 'e', 'g'],
        ['', '', '', 'b', '', 'g', ''],
        ['f', '', '', 'e', 'g', '', ''],
        ['a', '', '', '', '', 'd', ''],
        ['', '', 'g', '', '', '', 'a'],
        ['', 'g', '', 'd', '', '', ''],
        ['', 'c', 'd', '', 'a', 'b', '']
    ]
    
    for i in range(7):
        for j in range(7):
            if initial[i][j] and initial[i][j] != grid[i][j]:
                return False
    
    return True

def generate_solution():
    # Try each possible letter for minor diagonal
    for minor_letter in 'abcdefg':
        # Create base grid with minor diagonal filled
        grid = [['' for _ in range(7)] for _ in range(7)]
        
        # Fill minor diagonal
        for i in range(7):
            grid[i][6-i] = minor_letter
        
        # Fill known positions
        initial = [
            ['', 'd', '', '', 'b', 'e', 'g'],
            ['', '', '', 'b', '', 'g', ''],
            ['f', '', '', 'e', 'g', '', ''],
            ['a', '', '', '', '', 'd', ''],
            ['', '', 'g', '', '', '', 'a'],
            ['', 'g', '', 'd', '', '', ''],
            ['', 'c', 'd', '', 'a', 'b', '']
        ]
        
        valid = True
        for i in range(7):
            for j in range(7):
                if initial[i][j]:
                    if j == 6-i:  # On minor diagonal
                        if initial[i][j] != minor_letter:
                            valid = False
                            break
                    else:
                        grid[i][j] = initial[i][j]
            if not valid:
                break
        
        if not valid:
            continue
        
        # Try to fill remaining positions systematically
        remaining_positions = []
        for i in range(7):
            for j in range(7):
                if grid[i][j] == '' and j != 6-i:
                    remaining_positions.append((i, j))
        
        # For each position, try each possible letter
        def fill_remaining(pos_idx):
            if pos_idx == len(remaining_positions):
                return is_valid_complete_solution(grid)
            
            i, j = remaining_positions[pos_idx]
            used_in_row = set(grid[i])
            used_in_col = set(grid[r][j] for r in range(7))
            
            for letter in 'abcdefg':
                if letter not in used_in_row and letter not in used_in_col:
                    grid[i][j] = letter
                    if fill_remaining(pos_idx + 1):
                        return True
                    grid[i][j] = ''
            return False
        
        if fill_remaining(0):
            return grid
    
    return None

solution = generate_solution()
if solution:
    print_solution(solution)
else:
    print("No solution found")