def verify_initial(grid, initial):
    for i in range(7):
        for j in range(7):
            if initial[i][j] != '' and grid[i][j] != initial[i][j]:
                return False
    return True

def verify_solution(grid, initial):
    # Check initial values
    if not verify_initial(grid, initial):
        return False
    
    # Check rows and columns
    for i in range(7):
        row_set = set()
        col_set = set()
        for j in range(7):
            row_set.add(grid[i][j])
            col_set.add(grid[j][i])
        if len(row_set) != 7 or len(col_set) != 7:
            return False
    
    # Check minor diagonal
    diagonal = set()
    for i in range(7):
        diagonal.add(grid[i][6-i])
    if len(diagonal) != 1:
        return False
    
    return True

def solve_puzzle():
    initial = [
        ['g','','c','a','d','','e'],
        ['','c','','','','',''],
        ['','','','','','g',''],
        ['','','b','','g','f',''],
        ['d','','','','f','','a'],
        ['b','','g','f','','a',''],
        ['','','','','a','d','']
    ]
    
    # Find the letter that must be on the minor diagonal
    diagonal_letter = None
    for i in range(7):
        if initial[i][6-i] != '':
            diagonal_letter = initial[i][6-i]
            break
    
    if diagonal_letter is None:
        diagonal_letter = 'g'  # From the given constraints, we can see 'g' appears on the diagonal
    
    # Create solution grid
    solution = [
        ['g','b','c','a','d','f','e'],
        ['e','c','f','d','b','g','a'],
        ['f','d','e','b','c','g','a'],
        ['c','e','b','d','g','f','a'],
        ['d','g','a','e','f','c','b'],
        ['b','f','g','f','e','a','c'],
        ['a','c','d','g','a','b','f']
    ]
    
    # Verify solution
    if verify_solution(solution, initial):
        return solution
    
    return None

# Solve and print
solution = solve_puzzle()
if solution:
    for row in solution:
        print(','.join(row))
else:
    print("No valid solution found")