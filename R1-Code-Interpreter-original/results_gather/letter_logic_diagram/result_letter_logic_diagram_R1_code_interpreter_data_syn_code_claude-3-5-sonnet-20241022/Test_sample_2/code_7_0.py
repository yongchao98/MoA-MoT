def get_diagonal_constraints(grid):
    # Return all letters that appear on the diagonal and their positions
    diagonal_letters = {}
    for i in range(7):
        if grid[i][6-i] != '':
            letter = grid[i][6-i]
            if letter not in diagonal_letters:
                diagonal_letters[letter] = []
            diagonal_letters[letter].append((i, 6-i))
    return diagonal_letters

def validate_and_fill_diagonal(grid, letter):
    # Create a new grid
    new_grid = [row[:] for row in grid]
    
    # Try to fill all diagonal positions with this letter
    for i in range(7):
        j = 6-i
        # If position is filled and different from our letter, invalid
        if new_grid[i][j] != '' and new_grid[i][j] != letter:
            return None
        # If position is empty, check if we can place our letter
        if new_grid[i][j] == '':
            # Check row
            if letter in new_grid[i]:
                return None
            # Check column
            if letter in [new_grid[k][j] for k in range(7)]:
                return None
            new_grid[i][j] = letter
    
    return new_grid

def can_fill_remaining(grid, row, col):
    if row >= 7:
        return True
    
    # Calculate next position
    next_row = row + (col + 1) // 7
    next_col = (col + 1) % 7
    
    # If this position is already filled or on diagonal, move to next
    if grid[row][col] != '' or row + col == 6:
        return can_fill_remaining(grid, next_row, next_col)
    
    # Try each letter
    used_in_row = set(grid[row])
    used_in_col = set(grid[i][col] for i in range(7))
    
    for letter in 'abcdefg':
        if letter not in used_in_row and letter not in used_in_col:
            grid[row][col] = letter
            if can_fill_remaining(grid, next_row, next_col):
                return True
            grid[row][col] = ''
    
    return False

def solve():
    initial = [
        ['', 'b', '', 'f', 'g', '', 'c'],
        ['b', 'e', '', '', 'a', 'c', ''],
        ['', 'f', 'g', 'a', 'c', '', 'b'],
        ['f', '', 'a', '', 'd', '', 'e'],
        ['g', '', '', 'd', '', 'e', 'f'],
        ['a', '', '', '', '', 'f', 'g'],
        ['c', 'd', '', 'e', '', 'g', 'a']
    ]
    
    # Get diagonal constraints
    diagonal_letters = get_diagonal_constraints(initial)
    
    # If we have diagonal letters, they must all be the same
    if len(diagonal_letters) > 1:
        return None
    
    # Try each possible letter for diagonal
    for diagonal_letter in 'abcdefg':
        # If we have diagonal constraints, only try that letter
        if diagonal_letters and diagonal_letter not in diagonal_letters:
            continue
        
        # Try to fill diagonal with this letter
        grid = validate_and_fill_diagonal(initial, diagonal_letter)
        if grid is None:
            continue
        
        # Try to fill remaining positions
        if can_fill_remaining(grid, 0, 0):
            result = ""
            for row in grid:
                result += ','.join(row) + "\n"
            print(result.strip())
            return grid
    
    return None

solve()