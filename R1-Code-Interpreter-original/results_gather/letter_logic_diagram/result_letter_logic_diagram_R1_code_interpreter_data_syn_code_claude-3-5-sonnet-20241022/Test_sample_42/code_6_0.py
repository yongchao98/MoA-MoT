original = [
    ['', 'b', '', '', 'g', 'd', 'c'],
    ['', 'a', 'f', 'g', '', 'c', 'e'],
    ['a', 'f', '', '', 'c', 'e', 'b'],
    ['', 'g', 'd', '', 'e', '', ''],
    ['g', 'd', 'c', 'e', 'b', '', 'f'],
    ['', 'c', 'e', '', 'a', '', ''],
    ['', '', 'b', 'a', 'f', '', '']
]

solution = [
    ['e','b','a','f','g','d','c'],
    ['b','a','f','g','d','c','e'],
    ['a','f','g','d','c','e','b'],
    ['f','g','d','c','e','b','a'],
    ['g','d','c','e','b','a','f'],
    ['d','c','e','b','a','f','g'],
    ['c','e','b','a','f','g','d']
]

# Verify all constraints
def verify_all():
    # 1. Check pre-filled cells
    for i in range(7):
        for j in range(7):
            if original[i][j] != '' and original[i][j] != solution[i][j]:
                return False, f"Pre-filled cell mismatch at ({i},{j})"
    
    # 2. Check rows
    for row in solution:
        if sorted(row) != ['a','b','c','d','e','f','g']:
            return False, "Row constraint violated"
    
    # 3. Check columns
    for j in range(7):
        col = [solution[i][j] for i in range(7)]
        if sorted(col) != ['a','b','c','d','e','f','g']:
            return False, "Column constraint violated"
    
    # 4. Check minor diagonal
    diagonal = [solution[i][6-i] for i in range(7)]
    if not all(x == 'c' for x in diagonal):
        return False, "Minor diagonal constraint violated"
    
    return True, "All constraints satisfied"

result, message = verify_all()
print(f"Solution verification: {result}")
print(f"Message: {message}")