def fill_grid(grid):
    from collections import Counter
    
    # Determine the letters already used in each row and column
    rows = [set(row) for row in grid]
    cols = [set(col) for col in zip(*grid)]
    
    # Determine the letter for the minor diagonal
    possible_letters = set('abcdefg')
    for i, j in [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]:
        if grid[i][j] != '':
            possible_letters.intersection_update(grid[i][j])
    
    # Choose a letter for the minor diagonal
    minor_diagonal_letter = possible_letters.pop()
    
    # Fill the minor diagonal
    for i, j in [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]:
        grid[i][j] = minor_diagonal_letter
        rows[i].add(minor_diagonal_letter)
        cols[j].add(minor_diagonal_letter)
    
    # Fill the rest of the grid
    for i in range(7):
        for j in range(7):
            if grid[i][j] == '':
                for letter in 'abcdefg':
                    if letter not in rows[i] and letter not in cols[j]:
                        grid[i][j] = letter
                        rows[i].add(letter)
                        cols[j].add(letter)
                        break
    
    return grid

# Initial grid
grid = [
    ['b', '', 'e', '', '', '', ''],
    ['', '', '', '', 'f', '', ''],
    ['', '', '', 'f', 'g', '', ''],
    ['c', '', 'f', 'g', 'b', 'd', ''],
    ['', 'f', '', 'b', 'd', '', 'c'],
    ['', 'g', '', 'd', '', '', 'a'],
    ['g', '', '', '', '', '', 'f']
]

# Fill the grid
filled_grid = fill_grid(grid)

# Print the filled grid
for row in filled_grid:
    print(','.join(row))