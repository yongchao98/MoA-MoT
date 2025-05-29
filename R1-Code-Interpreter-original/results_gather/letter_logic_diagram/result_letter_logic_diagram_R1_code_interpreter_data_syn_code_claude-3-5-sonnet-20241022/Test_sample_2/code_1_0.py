def print_grid(grid):
    for row in grid:
        print(','.join(row))

def get_row_constraints(grid):
    # For each row, get the set of letters already used
    constraints = []
    for row in grid:
        constraints.append(set(x for x in row if x != ''))
    return constraints

def get_col_constraints(grid):
    # For each column, get the set of letters already used
    constraints = [set() for _ in range(7)]
    for i in range(7):
        for j in range(7):
            if grid[i][j] != '':
                constraints[j].add(grid[i][j])
    return constraints

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
    
    # Force 'b' as the diagonal letter (we can try this first as it's less constrained)
    working_grid = [row[:] for row in grid]
    
    # Fill diagonal with 'b'
    for i in range(7):
        if working_grid[i][6-i] != '' and working_grid[i][6-i] != 'b':
            return None
        working_grid[i][6-i] = 'b'
    
    row_constraints = get_row_constraints(working_grid)
    col_constraints = get_col_constraints(working_grid)
    
    def can_place(letter, row, col):
        if letter in row_constraints[row]:
            return False
        if letter in col_constraints[col]:
            return False
        return True
    
    def find_empty():
        for i in range(7):
            for j in range(7):
                if working_grid[i][j] == '' and (i + j != 6):  # Skip diagonal
                    return (i, j)
        return None
    
    def solve_remaining():
        pos = find_empty()
        if not pos:
            return True
        
        row, col = pos
        for letter in 'acdefg':  # Try all letters except 'b' (diagonal)
            if can_place(letter, row, col):
                working_grid[row][col] = letter
                row_constraints[row].add(letter)
                col_constraints[col].add(letter)
                
                if solve_remaining():
                    return True
                    
                working_grid[row][col] = ''
                row_constraints[row].remove(letter)
                col_constraints[col].remove(letter)
        
        return False
    
    if solve_remaining():
        result = ""
        for row in working_grid:
            result += ','.join(row) + "\n"
        print(result.strip())
        return True
    
    return False

solve()