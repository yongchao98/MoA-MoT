def analyze_grid(grid):
    # For each position on the diagonal, determine what letters are possible
    possible_letters = {i: set('abcdefg') for i in range(7)}
    
    # Process pre-filled cells
    for i in range(7):
        for j in range(7):
            if grid[i][j] != '':
                letter = grid[i][j]
                # Remove this letter as possibility from all positions in same row
                for col in range(7):
                    if col != j and col + i == 6:  # if it's a diagonal position
                        if letter in possible_letters[i]:
                            possible_letters[i].remove(letter)
                
                # Remove this letter as possibility from all positions in same column
                for row in range(7):
                    if row != i and row + j == 6:  # if it's a diagonal position
                        if letter in possible_letters[row]:
                            possible_letters[row].remove(letter)
                            
                # If this is a diagonal position, this is the only possible letter
                if i + j == 6:
                    possible_letters[i] = {letter}
    
    # Find common letters across all diagonal positions
    common_letters = set('abcdefg')
    for pos_set in possible_letters.values():
        common_letters &= pos_set
    
    return list(common_letters)[0] if common_letters else None

def create_solution_with_diagonal(grid, diagonal_letter):
    # Create new grid with diagonal letter filled
    new_grid = [row[:] for row in grid]
    for i in range(7):
        j = 6-i
        if new_grid[i][j] == '':
            new_grid[i][j] = diagonal_letter
        elif new_grid[i][j] != diagonal_letter:
            return None
    
    # Fill remaining positions
    def solve(pos=0):
        if pos == 49:  # 7x7 = 49 positions
            return True
            
        row, col = pos // 7, pos % 7
        
        # Skip filled positions and diagonal positions
        if new_grid[row][col] != '' or row + col == 6:
            return solve(pos + 1)
        
        # Try each letter
        for letter in 'abcdefg':
            if letter != diagonal_letter:  # Don't use diagonal letter elsewhere
                # Check if letter can be placed
                valid = True
                # Check row
                for j in range(7):
                    if new_grid[row][j] == letter:
                        valid = False
                        break
                # Check column
                if valid:
                    for i in range(7):
                        if new_grid[i][col] == letter:
                            valid = False
                            break
                
                if valid:
                    new_grid[row][col] = letter
                    if solve(pos + 1):
                        return True
                    new_grid[row][col] = ''
        
        return False
    
    if solve():
        return new_grid
    return None

# Initial grid
grid = [
    ['a','e','f','g','','c','d'],
    ['e','f','g','b','','d','a'],
    ['','g','b','','','','e'],
    ['','b','','','a','','f'],
    ['','c','d','','e','','g'],
    ['c','d','a','e','f','g','b'],
    ['d','a','e','f','','b','']
]

# Find the required diagonal letter
diagonal_letter = analyze_grid(grid)
if diagonal_letter:
    solution = create_solution_with_diagonal(grid, diagonal_letter)
    if solution:
        for row in solution:
            print(','.join(row))
    else:
        print("No solution exists")
else:
    print("Could not determine diagonal letter")