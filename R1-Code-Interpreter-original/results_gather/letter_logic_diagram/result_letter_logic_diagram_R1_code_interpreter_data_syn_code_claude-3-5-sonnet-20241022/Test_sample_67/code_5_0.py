def print_solution(grid):
    print('<<<')
    for row in grid:
        print(','.join(row))
    print('>>>')

def is_valid_complete(grid):
    # Check rows and columns
    for i in range(7):
        row_set = set(grid[i])
        col_set = set(grid[j][i] for j in range(7))
        if len(row_set) != 7 or len(col_set) != 7:
            return False
    
    # Check diagonal
    diagonal = [grid[i][6-i] for i in range(7)]
    if len(set(diagonal)) != 1:
        return False
    
    # Check pre-filled positions
    initial = [
        ['', 'c', 'b', '', '', 'd', ''],
        ['c', '', 'g', 'e', 'd', 'f', 'a'],
        ['b', '', '', '', '', '', 'c'],
        ['', '', '', '', 'a', '', ''],
        ['', '', 'f', '', '', '', ''],
        ['', 'f', '', 'c', '', 'g', ''],
        ['f', '', '', 'b', 'g', '', '']
    ]
    
    for i in range(7):
        for j in range(7):
            if initial[i][j] and initial[i][j] != grid[i][j]:
                return False
    
    return True

def generate_valid_grid():
    # Try each possible letter for the diagonal
    for diagonal_letter in 'abcdefg':
        # Initialize grid with empty cells
        grid = [['' for _ in range(7)] for _ in range(7)]
        
        # Fill in the pre-filled positions
        prefilled = [
            ['', 'c', 'b', '', '', 'd', ''],
            ['c', '', 'g', 'e', 'd', 'f', 'a'],
            ['b', '', '', '', '', '', 'c'],
            ['', '', '', '', 'a', '', ''],
            ['', '', 'f', '', '', '', ''],
            ['', 'f', '', 'c', '', 'g', ''],
            ['f', '', '', 'b', 'g', '', '']
        ]
        
        # Check if diagonal_letter conflicts with pre-filled positions
        can_use = True
        for i in range(7):
            j = 6-i
            if prefilled[i][j] and prefilled[i][j] != diagonal_letter:
                can_use = False
                break
        
        if not can_use:
            continue
            
        # Fill in pre-filled positions and diagonal
        for i in range(7):
            for j in range(7):
                if prefilled[i][j]:
                    grid[i][j] = prefilled[i][j]
                elif i + j == 6:  # diagonal position
                    grid[i][j] = diagonal_letter
        
        # Try to fill remaining positions
        letters = set('abcdefg')
        
        # Fill row by row
        for i in range(7):
            available = letters.copy()
            # Remove used letters in this row
            for j in range(7):
                if grid[i][j]:
                    available.discard(grid[i][j])
            
            # Fill empty positions
            for j in range(7):
                if not grid[i][j]:
                    # Try each available letter
                    for letter in list(available):
                        # Check if letter can be used in this column
                        used_in_col = set(grid[k][j] for k in range(7) if grid[k][j])
                        if letter not in used_in_col:
                            grid[i][j] = letter
                            available.remove(letter)
                            break
                    
                    if not grid[i][j]:  # If we couldn't fill this position
                        break
            
            if '' in grid[i]:  # If we couldn't fill this row
                break
        
        # Check if we found a valid solution
        if is_valid_complete(grid):
            return grid
    
    return None

# Generate and print solution
solution = generate_valid_grid()
if solution:
    print_solution(solution)
else:
    print("No solution exists")