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

# Let's verify if the solution matches pre-filled cells
def verify_prefilled():
    for i in range(7):
        for j in range(7):
            if original[i][j] != ' ' and original[i][j] != solution[i][j]:
                print(f"Mismatch at position ({i},{j})")
                return False
    return True

print("Solution matches pre-filled cells:", verify_prefilled())