def analyze_puzzle():
    initial_grid = [
        ['', '', 'd', 'g', '', 'b', 'f'],
        ['', 'd', '', '', '', 'f', ''],
        ['d', '', 'g', '', '', 'f', ''],
        ['g', '', 'b', 'f', 'e', 'c', 'd'],
        ['a', 'b', 'f', '', 'c', 'd', ''],
        ['', '', '', 'c', 'd', 'g', ''],
        ['f', '', 'c', 'd', 'g', 'a', 'b']
    ]
    
    # Check diagonal cells
    diagonal_letters = set()
    for i in range(7):
        j = 6 - i
        if initial_grid[i][j] != '':
            diagonal_letters.add(initial_grid[i][j])
    
    print("Letters on diagonal:", diagonal_letters)
    
    # Check each row's pre-filled letters
    for i in range(7):
        letters = [x for x in initial_grid[i] if x != '']
        if len(letters) != len(set(letters)):
            print(f"Row {i+1} has duplicates:", letters)
        print(f"Row {i+1} letters:", letters)
    
    # Check each column's pre-filled letters
    for j in range(7):
        letters = [initial_grid[i][j] for i in range(7) if initial_grid[i][j] != '']
        if len(letters) != len(set(letters)):
            print(f"Column {j+1} has duplicates:", letters)
        print(f"Column {j+1} letters:", letters)
    
    # Count occurrences of each letter
    letter_count = {}
    for letter in 'abcdefg':
        count = sum(row.count(letter) for row in initial_grid)
        letter_count[letter] = count
    print("\nLetter counts:", letter_count)
    
    # Check diagonal constraints
    diagonal_positions = [(i, 6-i) for i in range(7)]
    print("\nDiagonal positions and their values:")
    for i, j in diagonal_positions:
        print(f"Position ({i},{j}): {initial_grid[i][j]}")

analyze_puzzle()