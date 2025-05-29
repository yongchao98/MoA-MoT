def verify_solution(grid, original):
    # Check pre-filled cells
    for i in range(7):
        for j in range(7):
            if original[i][j] and original[i][j] != grid[i][j]:
                return False
    
    # Check rows and columns
    letters = set('abcdefg')
    for i in range(7):
        if set(grid[i]) != letters:  # Check row
            return False
        if set(grid[j][i] for j in range(7)) != letters:  # Check column
            return False
    
    # Check minor diagonal
    diag = [grid[i][6-i] for i in range(7)]
    if len(set(diag)) != 1:  # All must be same letter
        return False
    
    return True

def find_diagonal_letter(original):
    # Count frequency of letters on minor diagonal in original
    freq = {}
    for i in range(7):
        if original[i][6-i]:
            letter = original[i][6-i]
            freq[letter] = freq.get(letter, 0) + 1
    
    # Return letter with highest frequency
    if freq:
        return max(freq.items(), key=lambda x: x[1])[0]
    return None

def solve_puzzle():
    original = [
        ['f','','','c','a','','g'],
        ['e','b','c','a','d','g',''],
        ['','c','','','g','f','e'],
        ['c','a','','','','e','b'],
        ['a','d','g','','e','','c'],
        ['d','g','','e','','c',''],
        ['','','','b','','a','']
    ]
    
    # Create empty grid
    grid = [[''] * 7 for _ in range(7)]
    
    # Copy pre-filled cells
    for i in range(7):
        for j in range(7):
            grid[i][j] = original[i][j]
    
    # Find and fill diagonal letter
    diag_letter = find_diagonal_letter(original)
    if not diag_letter:
        return None
    
    # Fill all diagonal positions with diag_letter
    for i in range(7):
        if not grid[i][6-i]:
            grid[i][6-i] = diag_letter
    
    def can_place(row, col, letter):
        if grid[row][col]:  # Skip filled cells
            return False
        
        # Check row
        if letter in grid[row]:
            return False
        
        # Check column
        if letter in [grid[i][col] for i in range(7)]:
            return False
        
        return True
    
    def solve(pos=0):
        if pos == 49:
            return verify_solution(grid, original)
        
        row, col = pos // 7, pos % 7
        
        # Skip filled cells
        if grid[row][col]:
            return solve(pos + 1)
        
        for letter in 'abcdefg':
            if can_place(row, col, letter):
                grid[row][col] = letter
                if solve(pos + 1):
                    return True
                grid[row][col] = ''
        
        return False
    
    return grid if solve() else None

# Solve and print
solution = solve_puzzle()
if solution:
    for row in solution:
        print(','.join(row))
else:
    print("No solution found")