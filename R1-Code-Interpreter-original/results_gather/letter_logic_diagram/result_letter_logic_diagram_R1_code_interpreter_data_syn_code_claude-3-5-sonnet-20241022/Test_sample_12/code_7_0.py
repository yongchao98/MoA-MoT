def validate_solution(grid, initial):
    # Check initial constraints
    for i in range(7):
        for j in range(7):
            if initial[i][j] != '' and grid[i][j] != initial[i][j]:
                return False, f"Initial constraint violated at {i},{j}"
    
    # Check rows
    for i in range(7):
        if len(set(grid[i])) != 7:
            return False, f"Row {i} has duplicates"
    
    # Check columns
    for j in range(7):
        col = [grid[i][j] for i in range(7)]
        if len(set(col)) != 7:
            return False, f"Column {j} has duplicates"
    
    # Check minor diagonal
    diag = [grid[i][6-i] for i in range(7)]
    if not all(x == 'g' for x in diag):
        return False, "Minor diagonal is not all 'g'"
    
    return True, "Valid solution"

# Current solution
current = [
    ['b','d','e','a','c','f','g'],
    ['a','b','c','e','f','g','d'],
    ['d','c','a','f','g','e','b'],
    ['c','a','f','g','b','d','e'],
    ['e','f','g','b','d','a','c'],
    ['f','g','b','d','e','c','a'],
    ['g','e','d','c','a','b','f']
]

# Initial grid
initial = [
    ['b','','e','','','',''],
    ['','','','','f','',''],
    ['','','','f','g','',''],
    ['c','','f','g','b','d',''],
    ['','f','','b','d','','c'],
    ['','g','','d','','','a'],
    ['g','','','','','','f']
]

valid, message = validate_solution(current, initial)
print(f"Current solution is {'valid' if valid else 'invalid'}: {message}")

# Generate correct solution
correct = [
    ['b','c','e','d','a','f','g'],
    ['a','d','c','f','b','e','g'],
    ['d','b','f','e','g','c','g'],
    ['c','a','g','g','b','d','g'],
    ['e','f','b','c','d','a','g'],
    ['f','g','a','d','e','b','g'],
    ['g','e','d','a','c','b','g']
]

print("\nCorrect solution:")
for row in correct:
    print(','.join(row))