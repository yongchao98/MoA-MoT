original = [
    [' ','b',' ',' ','g','d','c'],
    [' ','a','f','g',' ','c','e'],
    ['a','f',' ',' ','c','e','b'],
    [' ','g','d',' ','e',' ',' '],
    ['g','d','c','e','b',' ','f'],
    [' ','c','e',' ','a',' ',' '],
    [' ',' ','b','a','f',' ',' ']
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

def verify_solution():
    # Check pre-filled cells
    for i in range(7):
        for j in range(7):
            if original[i][j] != ' ' and original[i][j] != solution[i][j]:
                print(f"Mismatch at ({i},{j}): expected {original[i][j]}, got {solution[i][j]}")
                return False
    return True

print("Solution matches pre-filled cells:", verify_solution())