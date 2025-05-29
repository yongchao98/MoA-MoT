def print_solution(grid):
    for row in grid:
        print(','.join(row))

def check_conflicts(grid, row, col, letter):
    # Check row
    if letter in grid[row]:
        return True
    
    # Check column
    if letter in [grid[i][col] for i in range(7)]:
        return True
    
    return False

def solve():
    # Initial grid
    grid = [
        ['', 'b', '', 'f', 'g', '', 'c'],
        ['b', 'e', '', '', 'a', 'c', ''],
        ['', 'f', 'g', 'a', 'c', '', 'b'],
        ['f', '', 'a', '', 'd', '', 'e'],
        ['g', '', '', 'd', '', 'e', 'f'],
        ['a', '', '', '', '', 'f', 'g'],
        ['c', 'd', '', 'e', '', 'g', 'a']
    ]
    
    # First, find what letters are already on the minor diagonal
    diagonal_letters = set()
    for i in range(7):
        if grid[i][6-i] != '':
            diagonal_letters.add(grid[i][6-i])
    
    # Try each possible letter for the diagonal
    for diagonal_letter in 'abcdefg':
        # Skip if this letter conflicts with pre-filled diagonal positions
        if diagonal_letters and diagonal_letter not in diagonal_letters:
            continue
            
        # Create a fresh copy of the grid
        test_grid = [row[:] for row in grid]
        
        # Fill all diagonal positions with this letter
        valid = True
        for i in range(7):
            if test_grid[i][6-i] == '':
                # Check if this placement creates conflicts
                if check_conflicts(test_grid, i, 6-i, diagonal_letter):
                    valid = False
                    break
                test_grid[i][6-i] = diagonal_letter
            elif test_grid[i][6-i] != diagonal_letter:
                valid = False
                break
                
        if not valid:
            continue
        
        def solve_remaining(pos=0):
            if pos >= 49:
                return True
            
            row = pos // 7
            col = pos % 7
            
            # Skip filled positions and diagonal
            while pos < 49 and (test_grid[row][col] != '' or row + col == 6):
                pos += 1
                if pos < 49:
                    row = pos // 7
                    col = pos % 7
            
            if pos >= 49:
                return True
            
            # Try each letter except diagonal_letter
            for letter in [l for l in 'abcdefg' if l != diagonal_letter]:
                if not check_conflicts(test_grid, row, col, letter):
                    test_grid[row][col] = letter
                    if solve_remaining(pos + 1):
                        return True
                    test_grid[row][col] = ''
            
            return False
        
        if solve_remaining():
            # Verify solution
            for i in range(7):
                if len(set(test_grid[i])) != 7:  # Check rows
                    continue
                if len(set(test_grid[j][i] for j in range(7))) != 7:  # Check columns
                    continue
            
            result = ""
            for row in test_grid:
                result += ','.join(row) + "\n"
            print(result.strip())
            return True
    
    return False

solve()