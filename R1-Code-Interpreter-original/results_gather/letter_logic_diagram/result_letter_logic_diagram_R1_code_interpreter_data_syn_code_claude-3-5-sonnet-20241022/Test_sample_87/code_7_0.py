def print_solution(grid):
    for row in grid:
        print(','.join(row))

def is_valid_solution(grid):
    # Check rows and columns
    for i in range(7):
        row_set = set(grid[i])
        col_set = set(grid[j][i] for j in range(7))
        if len(row_set) != 7 or len(col_set) != 7:
            return False
    
    # Check diagonal (must be all same letter)
    diagonal = [grid[i][6-i] for i in range(7)]
    if len(set(diagonal)) != 1:
        return False
    
    # Check pre-filled values
    initial = [
        ['e', 'a', 'c', '', '', '', ''],
        ['', '', 'g', 'f', '', 'd', 'e'],
        ['c', 'g', 'f', 'b', 'd', '', ''],
        ['', '', 'b', 'd', 'e', '', 'c'],
        ['', '', '', '', '', '', ''],
        ['', 'd', 'e', '', 'c', 'g', 'f'],
        ['d', '', '', '', 'g', 'f', '']
    ]
    
    for i in range(7):
        for j in range(7):
            if initial[i][j] != '' and grid[i][j] != initial[i][j]:
                return False
    
    return True

def solve(grid, diagonal_letter):
    def find_next_empty():
        # First fill diagonal positions
        for i in range(7):
            if i + (6-i) == 6 and grid[i][6-i] == '':
                return (i, 6-i)
        # Then other positions
        for i in range(7):
            for j in range(7):
                if grid[i][j] == '':
                    return (i, j)
        return None

    def is_valid_placement(row, col, letter):
        # Check row
        for j in range(7):
            if grid[row][j] == letter:
                return False
        
        # Check column
        for i in range(7):
            if grid[i][col] == letter:
                return False
        
        # Check diagonal requirement
        if row + col == 6 and letter != diagonal_letter:
            return False
        
        return True

    pos = find_next_empty()
    if not pos:
        return is_valid_solution(grid)
    
    row, col = pos
    letters = [diagonal_letter] if row + col == 6 else [c for c in 'abcdefg']
    
    for letter in letters:
        if is_valid_placement(row, col, letter):
            grid[row][col] = letter
            if solve(grid, diagonal_letter):
                return True
            grid[row][col] = ''
    
    return False

# Initial grid
grid = [
    ['e', 'a', 'c', '', '', '', ''],
    ['', '', 'g', 'f', '', 'd', 'e'],
    ['c', 'g', 'f', 'b', 'd', '', ''],
    ['', '', 'b', 'd', 'e', '', 'c'],
    ['', '', '', '', '', '', ''],
    ['', 'd', 'e', '', 'c', 'g', 'f'],
    ['d', '', '', '', 'g', 'f', '']
]

print('<<<')
found_solution = False
for diagonal_letter in 'abcdefg':
    grid_copy = [row[:] for row in grid]
    if solve(grid_copy, diagonal_letter):
        print_solution(grid_copy)
        found_solution = True
        break

if not found_solution:
    print("No solution exists")
print('>>>')